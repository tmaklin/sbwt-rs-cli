use super::kmer::*;
use super::kmer_chunk::KmerChunk;
use crate::tempfile::{TempFile, TempFileManager};
use crate::util::DNA_ALPHABET;
use std::io::{BufWriter, Seek, Write};
use std::borrow::BorrowMut;

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

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

pub fn split_to_bins<const B: usize, IN: crate::SeqStream + Send>(mut seqs: IN, k: usize, _mem_gb: usize, _n_threads: usize, _dedup_batches: bool, temp_file_manager: &mut TempFileManager) -> Vec<TempFile>{

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

    let mut buf = Vec::<Box<[u8]>>::new();
    while let Some(seq) = seqs.stream_next() {
        let mut seq_copy = seq.to_owned();
        seq_copy.reverse(); // Reverse to get colex sorting
        buf.push(seq_copy.into_boxed_slice());
    }

    // Create writers
    let mut bin_writers = colex_sorted_binmers(bin_prefix_len).into_iter().map(|binmer| {
        let name_prefix = format!("sbwt-temp-{}-", String::from_utf8(binmer).unwrap());
        temp_file_manager.create_new_file(&name_prefix, 8, ".bin")
    }).collect::<Vec<TempFile>>();

    let kmers: Vec<LongKmer::<B>> = buf.par_iter().map(|seq| {
        seq.windows(k).filter_map(|kmer| {
            LongKmer::<B>::from_ascii(kmer).ok()
        }).collect::<Vec<LongKmer::<B>>>()
    }).flatten().collect();

    kmers.iter().for_each(|kmer| {
        let bin_id = kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize; // Interpret nucleotides in base-4
        let bin_file = &mut bin_writers[bin_id];
        kmer.serialize(bin_file).unwrap(); // Todo: write all at once
    });

    // Return the TempFiles
    for w in bin_writers.iter_mut(){
        match w.file.rewind() {
            Ok(_) => (),
            Err(e) => panic!("Couldn't rewind file: {}", e)
        }
    }
    bin_writers

}

// Overwrite the files with sorted and deduplicates files. Returns back the files after overwriting.
pub fn par_sort_and_dedup_bin_files<const B: usize>(bin_files: &mut Vec<TempFile>, _mem_gb: usize, _n_threads: usize) {

    let mut files_and_sizes = bin_files.par_iter_mut().map(|f| {
        f
    }).collect::<Vec<&mut TempFile>>();
    files_and_sizes.sort_by_key(|f| f.avail_in());

    log::info!("Sorting k-mer bins");

    files_and_sizes.par_iter_mut().for_each(|f| {
        // Using debug log level as a more verbose info level
        let mut reader = std::io::BufReader::new(f.file.borrow_mut());
        let chunk = KmerChunk::<B>::load(&mut reader).unwrap();

        let mut chunk = chunk.sort();
        chunk.dedup();

        // Overwrite the file and seek to start
        f.file = std::io::Cursor::new(Vec::new());
        chunk.serialize(f.file.borrow_mut()).unwrap();
        f.file.seek(std::io::SeekFrom::Start(0)).unwrap();
    });
}

// The original files are deleted
pub fn concat_files(infiles: Vec<TempFile>, out_writer: &mut impl std::io::Write){
    let mut bw = BufWriter::new(out_writer);
    for mut fp in infiles {
        let mut reader = std::io::BufReader::new(&mut fp.file);
        std::io::copy(&mut reader, &mut bw).unwrap();
        drop(fp);
    }
    bw.flush().unwrap();
}

mod tests {
    #[test]
    fn test_colex_sorted_binmers(){
        let binmers = super::colex_sorted_binmers(2);
        let ans = vec![b"AA", b"CA", b"GA", b"TA", b"AC", b"CC", b"GC", b"TC", b"AG", b"CG", b"GG", b"TG", b"AT", b"CT", b"GT", b"TT"];
        assert_eq!(binmers, ans);
    }
}
