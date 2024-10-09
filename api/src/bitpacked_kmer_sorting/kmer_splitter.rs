use super::kmer::*;
use crate::tempfile::{TempFile, TempFileManager};
use crate::util::DNA_ALPHABET;
use std::io::{BufWriter, Write};

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSliceMut;
use rayon::prelude::ParallelSlice;

#[allow(dead_code)]
fn colex_sorted_binmers(bin_prefix_len: usize) -> Vec<Vec<u8>> {
    let mut binmers = Vec::<Vec<u8>>::new();
    for i in 0..(4_usize.pow(bin_prefix_len as u32)){
        let mut binmer = Vec::<u8>::new();
        let mut j = i;
        for _ in 0..bin_prefix_len{
            binmer.push(DNA_ALPHABET[j % 4]);
            j /= 4;
        }
        binmers.push(binmer);
    }
    binmers
}

pub fn split_to_bins<const B: usize, IN: crate::SeqStream + Send>(
    mut seqs: IN,
    k: usize,
    _mem_gb: usize,
    _n_threads: usize,
    _dedup_batches: bool,
    _temp_file_manager: &mut TempFileManager
) -> Vec<std::io::Cursor::<Vec<LongKmer::<B>>>> {

    // Suppose we have a memory budget of m bytes and t threads.
    // Suppose each k-mer takes s bytes and there are 64 bins.
    // Let b be the number of k-mers in each splitter thread bin buffer.
    // A splitter thread uses 64bs bytes
    // In total the splitter threads use 64bst threads.
    // So, we need:
    // 64bbt = m
    // b = m / (64bt)

    // Wrap to scope to be able to borrow seqs for the producer thread even when it's not 'static.
    let bin_prefix_len = 3_usize; // If you update this you must update all the logic below
    let n_bins = (4_usize).pow(bin_prefix_len as u32); // 64

    log::info!("Splitting k-mers into {} bins", n_bins);

    // This function has some weird implicit assumptions that break something if
    // the operations are done in another order??

    let mut buf = Vec::<Box<[u8]>>::new();
    while let Some(seq) = seqs.stream_next() {
        let mut seq_copy = seq.to_owned();
        seq_copy.reverse(); // Reverse to get colex sorting
        buf.push(seq_copy.into_boxed_slice());
    }

    let mut kmers = buf.par_iter().map(|seq| {
        seq.windows(k).filter_map(|bytes| {
            LongKmer::<B>::from_ascii(bytes).ok()
        }).collect::<Vec<LongKmer::<B>>>()
    }).flatten().collect::<Vec<LongKmer::<B>>>();

    kmers.par_sort_unstable_by_key(|kmer| {
        kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize
    });

    kmers.par_chunk_by(|&a, &b| {
        let x = a.get_from_left(0) as usize * 16 + a.get_from_left(1) as usize * 4 + a.get_from_left(2) as usize;
        let y = b.get_from_left(0) as usize * 16 + b.get_from_left(1) as usize * 4 + b.get_from_left(2) as usize;
        x == y}).map(|chunk| {
        std::io::Cursor::<Vec<LongKmer::<B>>>::new(chunk.to_vec())
    }).collect::<Vec<std::io::Cursor::<Vec<LongKmer::<B>>>>>()
}

// Overwrite the files with sorted and deduplicates files. Returns back the files after overwriting.
pub fn par_sort_and_dedup_bin_files<const B: usize>(bin_files: &mut Vec<std::io::Cursor::<Vec<LongKmer::<B>>>>, _mem_gb: usize, _n_threads: usize) {
    log::info!("Sorting k-mer bins");
    bin_files.par_iter_mut().for_each(|f| {
        f.get_mut().sort_unstable();
        f.get_mut().dedup();
        f.set_position(0);
    });
}

// The original files are deleted
#[allow(dead_code)]
pub fn concat_files(infiles: Vec<TempFile>, out_writer: &mut impl std::io::Write){
    let mut bw = BufWriter::new(out_writer);
    for mut fp in infiles {
        let mut reader = std::io::BufReader::new(&mut fp.file);
        std::io::copy(&mut reader, &mut bw).unwrap();
        drop(fp);
    }
    bw.flush().unwrap();
}

// The original files are deleted
pub fn concat_files_take<const B: usize>(infiles: &mut Vec<std::io::Cursor::<Vec<LongKmer::<B>>>>) -> TempFile {
    let mut writer: std::io::Cursor::<Vec<u8>> = std::io::Cursor::new(Vec::new());
    infiles.iter().for_each(|file| {
        file.get_ref().iter().for_each(|kmer| {
            kmer.serialize(&mut writer).expect("Serialized kmer to TempFile");
        })
    });
    writer.set_position(0);
    TempFile{ file: writer }
}

mod tests {
    #[test]
    fn test_colex_sorted_binmers(){
        let binmers = super::colex_sorted_binmers(2);
        let ans = vec![b"AA", b"CA", b"GA", b"TA", b"AC", b"CC", b"GC", b"TC", b"AG", b"CG", b"GG", b"TG", b"AT", b"CT", b"GT", b"TT"];
        assert_eq!(binmers, ans);
    }
}
