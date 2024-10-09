//! Build an [SbwtIndex] using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;
mod kmer;

use crate::{sbwt::{PrefixLookupTable, SbwtIndex}, streaming_index::LcsArray, subsetseq::SubsetSeq, tempfile::TempFileManager, util::DNA_ALPHABET};
use crate::bitpacked_kmer_sorting::kmer_splitter::concat_files_take;

/// Build using bitpacked k-mer sorting. See [SbwtIndexBuilder](crate::builder::SbwtIndexBuilder) for a wrapper with a more 
/// user-friendly interface. B is the number u64 words in a k-mer.
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, build_lcs: bool, temp_file_manager: &mut TempFileManager) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let sigma = DNA_ALPHABET.len();

    log::info!("Splitting k-mers into bins");
    let mut bin_files = kmer_splitter::split_to_bins::<B, IN>(seqs, k, mem_gb, n_threads, dedup_batches, temp_file_manager);

    log::info!("Sorting and deduplicating bins");
    kmer_splitter::par_sort_and_dedup_bin_files::<B>(&mut bin_files, mem_gb, n_threads);

    let mut kmers_file = concat_files_take(&mut bin_files);

    let n_kmers = kmers_file.get_ref().len();

    log::info!("{} distinct k-mers found", n_kmers);
    let mut required_dummies = dummies::get_sorted_dummies::<B>(&mut kmers_file, sigma, k, temp_file_manager);

    let n = n_kmers + required_dummies.get_ref().len();

    log::info!("Constructing the sbwt subset sequence");

    let (rawrows, lcs) = cursors::build_sbwt_bit_vectors::<B>(&mut kmers_file, &mut required_dummies, n, k, sigma, build_lcs);

    // Create the C array
    #[allow(non_snake_case)] // C-array is an established convention in BWT indexes
    let C: Vec<usize> = crate::util::get_C_array(&rawrows);

    log::info!("Building the subset rank structure");
    let mut subsetseq = SS::new_from_bit_vectors(rawrows.into_iter().map(simple_sds_sbwt::bit_vector::BitVector::from).collect());
    subsetseq.build_rank();
    let n_sets = subsetseq.len();
    let (mut index, lcs) = (SbwtIndex::<SS>::from_components(
        subsetseq,
        n_kmers,
        k,
        C,
        PrefixLookupTable::new_empty(n_sets))
                            , lcs.map(LcsArray::new));

    let lut = PrefixLookupTable::new(&index, 8);
    index.set_lookup_table(lut);
    (index, lcs)

}
