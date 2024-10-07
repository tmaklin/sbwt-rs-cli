use std::borrow::BorrowMut;
use std::io::{BufWriter, Write};

use super::kmer::LongKmer;
use crate::tempfile::TempFileManager;
use crate::tempfile::TempFile;
use simple_sds_sbwt::ops::SelectZero;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::raw_vector::*;
use rayon::prelude::*;

use crate::bitpacked_kmer_sorting::cursors::reset_reader_position;

#[allow(dead_code)]
struct NullReader{}

impl std::io::Read for NullReader{
    fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize>{
        Ok(0) // EoF
    }
}

// We take in a path and not a file object because we need multiple readers to the same file
pub fn get_sorted_dummies<const B: usize>(sorted_kmers: &mut TempFile, sigma: usize, k: usize, temp_file_manager: &mut TempFileManager) -> Vec<(LongKmer::<B>, u8)>{

    // Todo: I'm using dummy merger cursors with an empty dummy file. Should refactor things to manage without the
    // empty dummy file.

    // Number of k-mers in file
    let n = sorted_kmers.avail_in() as usize / LongKmer::<B>::byte_size();

    let mut has_predecessor = simple_sds_sbwt::raw_vector::RawVector::new();
    has_predecessor.resize(n, false);

    let mut emptyfile = temp_file_manager.create_new_file("empty-", 10, ".bin");
    let char_cursors = crate::bitpacked_kmer_sorting::cursors::init_char_cursor_positions::<B>(&mut emptyfile, sorted_kmers, k, sigma);

    let mut global_cursor = crate::bitpacked_kmer_sorting::cursors::DummyNodeMerger::new(
        &mut emptyfile,
        sorted_kmers,
        k,
    );

    let xs: Vec<(LongKmer::<B>, u8)> = global_cursor.borrow_mut().map(|x| {
        x
    }).collect();


    for c in 0..(sigma as u8) {
        reset_reader_position(&mut global_cursor,
                              char_cursors[c as usize].0.0,
                              char_cursors[c as usize].1.0,
                              char_cursors[c as usize].0.1 as usize,
                              char_cursors[c as usize].1.1 as usize);

        xs.iter().for_each(|(x, _)| {
            let xc = x.set_from_left(k-1, 0).right_shift(1).set_from_left(0, c);

            while global_cursor.peek().is_some(){
                match global_cursor.peek().unwrap().0.cmp(&xc) {
                    std::cmp::Ordering::Greater => {
                        break
                    },
                    std::cmp::Ordering::Equal => {
                        has_predecessor.set_bit(global_cursor.nondummy_position(), true);
                        global_cursor.next(); // Advance
                        break
                    },
                    std::cmp::Ordering::Less => {
                        global_cursor.next(); // Advance
                        // no break
                    }
                }
            }
        });
    }

    let iterable = BitVector::from(has_predecessor);
    let mut required_dummies: Vec::<(LongKmer::<B>, u8)> = iterable.zero_iter().par_bridge().map(|x| {
        let mut prefix = xs[x.1].0;
        (0..k).collect::<Vec<usize>>().iter().map(|i| {
            let len = k - i - 1;
            prefix = prefix.left_shift(1);
            (prefix, len as u8)
        }).collect::<Vec<(LongKmer::<B>, u8)>>()
    }).flatten().collect();

    // We always assume that the empty k-mer exists. This assumption is reflected in the C-arrya
    // later, which adds one "ghost dollar" count to all counts.
    required_dummies.push((LongKmer::<B>::from_ascii(b"").unwrap(), 0));

    required_dummies.par_sort();
    required_dummies.dedup();
    required_dummies.shrink_to_fit();

    required_dummies
    

}

pub fn write_to_disk<const B: usize>(dummies: Vec<(LongKmer::<B>, u8)>, writer: &mut TempFile){
    let mut bw = BufWriter::new(writer);
    for (kmer, len) in dummies.iter(){
        kmer.serialize(&mut bw).unwrap();
        bw.write_all(&[*len]).unwrap();
    }
    bw.flush().unwrap();
}
