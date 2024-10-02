use super::kmer::*;
use super::kmer_chunk::KmerChunk;
use crate::tempfile::{TempFile, TempFileManager};
use crate::util::DNA_ALPHABET;
use std::io::{BufWriter, Seek, Write};
use std::borrow::BorrowMut;

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

pub fn split_to_bins<const B: usize, IN: crate::SeqStream + Send>(mut seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, temp_file_manager: &mut TempFileManager) -> Vec<TempFile>{

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
    let mem_gb = (mem_gb as usize * (1_usize << 30) as usize);
    let encoder_bin_buf_size = mem_gb / ((n_bins as usize * LongKmer::<B>::byte_size() as usize) * n_threads as usize);

    log::info!("Splitting k-mers into {} bins", n_bins);
    let mut bin_buffers = vec![Vec::<LongKmer::<B>>::new(); n_bins];
    for buf in bin_buffers.iter_mut(){
        buf.reserve_exact(encoder_bin_buf_size.try_into().unwrap());
    }

    let mut buf = Vec::<Box<[u8]>>::new();
    while let Some(seq) = seqs.stream_next() {
	let mut seq_copy = seq.to_owned();
        seq_copy.reverse(); // Reverse to get colex sorting
	buf.push(seq_copy.into_boxed_slice());
    }

    for seq in &buf {
	for kmer in seq.windows(k){
            match LongKmer::<B>::from_ascii(kmer) {
		Ok(kmer) => {
                    let bin_id = kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize; // Interpret nucleotides in base-4
                    bin_buffers[bin_id].push(kmer);
		}
		Err(KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
		Err(KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
            }
	}
    }

    // Send remaining buffers
    if dedup_batches {
        for b in &mut bin_buffers{
            b.sort_unstable();
            b.dedup();
        }
    }

    // Create writers
    let mut bin_writers = Vec::<std::io::BufWriter::<TempFile>>::new();

    for binmer in colex_sorted_binmers(bin_prefix_len) {
        let name_prefix = format!("sbwt-temp-{}-", String::from_utf8(binmer).unwrap());
        let f = temp_file_manager.create_new_file(&name_prefix, 8, ".bin");
        bin_writers.push(BufWriter::new(f));
    }


    bin_buffers.iter().for_each(|batch| {
        if !batch.is_empty() {
            let bin_id = batch[0].get_from_left(0) as usize * 16 + batch[0].get_from_left(1) as usize * 4 + batch[0].get_from_left(2) as usize; // Intepret nucleotides in base-4
            let bin_file = &mut bin_writers[bin_id];
            for kmer in batch{
                kmer.serialize(bin_file).unwrap(); // Todo: write all at once
            }
        }
    });

    // Return the TempFiles
    let mut writers: Vec<TempFile> = bin_writers.into_iter().map(|w| w.into_inner().unwrap()).collect();
    for w in writers.iter_mut(){
        w.file.seek(std::io::SeekFrom::Start(0)).unwrap();
    }
    writers

}

// Overwrite the files with sorted and deduplicates files. Returns back the files after overwriting.
pub fn par_sort_and_dedup_bin_files<const B: usize>(bin_files: &mut Vec<TempFile>, mem_gb: usize, n_threads: usize) {

    let filesizes = bin_files.iter().map(|f| f.avail_in() as usize).collect::<Vec<usize>>();
    let mut files_and_sizes = bin_files.into_iter().enumerate().map(|(i, f)| (f, filesizes[i], i)).collect::<Vec<(&mut TempFile, usize, usize)>>();
        
    files_and_sizes.sort_by_key(|(_, size, _)| *size);

    let max_mem = mem_gb * (1_usize << 30);

    log::info!("Sorting k-mer bins");

    files_and_sizes.iter_mut().for_each(|(f, _s, _i)| {
        // Using debug log level as a more verbose info level
        let mut reader = std::io::BufReader::new(f.file.borrow_mut());
        let chunk = KmerChunk::<B>::load(&mut reader).unwrap();
        
        let mut chunk = chunk.sort();
        chunk.dedup();

        // Overwrite the file and seek to start
	f.file = std::io::Cursor::new(Vec::new());
        let chunk_out = std::io::BufWriter::new(f.file.borrow_mut());
        chunk.serialize(chunk_out).unwrap();
        f.file.seek(std::io::SeekFrom::Start(0)).unwrap();
    });
}

// The original files are deleted
pub fn concat_files(infiles: Vec<TempFile>, out_writer: &mut impl std::io::Write){
    let mut bw = BufWriter::new(out_writer);
    for mut fp in infiles {
        let mut reader = std::io::BufReader::new(&mut fp.file);
        std::io::copy(&mut reader, &mut bw).unwrap();
        // fp is dropped here, which deletes the file
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
