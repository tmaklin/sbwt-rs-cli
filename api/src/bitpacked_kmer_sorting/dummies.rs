use std::io::{BufWriter, Write};

use super::kmer::LongKmer;
use crate::tempfile::TempFileManager;
use crate::tempfile::TempFile;
use simple_sds_sbwt::ops::SelectZero;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::raw_vector::*;
use rayon::prelude::*;

use crate::bitpacked_kmer_sorting::cursors::DummyNodeMerger;
use crate::bitpacked_kmer_sorting::cursors::split_global_cursor;

#[allow(dead_code)]
struct NullReader{}

impl std::io::Read for NullReader{
    fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize>{
        Ok(0) // EoF
    }
}

pub fn get_set_bits<const B: usize>(
    kmers: &Vec<(LongKmer::<B>, u8)>,
    cursor: &mut DummyNodeMerger<std::io::Cursor::<Vec<u8>>, B>,
    k: usize,
    c: u8,
) -> Vec<usize> {
    let mut set_bits: Vec<usize> = Vec::new();
    kmers.iter().for_each(|(x, _)| {
        let xc = x.set_from_left(k-1, 0).right_shift(1).set_from_left(0, c as u8);

        while cursor.peek().is_some(){
            match cursor.peek().unwrap().0.cmp(&xc) {
                std::cmp::Ordering::Greater => {
                    break
                },
                std::cmp::Ordering::Equal => {
                    set_bits.push(cursor.nondummy_position());
                    // has_predecessor.set_bit(cursor.nondummy_position(), true);
                    cursor.next(); // Advance
                    break
                },
                std::cmp::Ordering::Less => {
                    cursor.next(); // Advance
                    // no break
                }
            }
        }
    });
    set_bits
}

// We take in a path and not a file object because we need multiple readers to the same file
pub fn get_sorted_dummies<const B: usize>(sorted_kmers: &mut TempFile, sigma: usize, k: usize, temp_file_manager: &mut TempFileManager) -> Vec<(LongKmer::<B>, u8)>{
    // Number of k-mers in file
    let n = sorted_kmers.avail_in() as usize / LongKmer::<B>::byte_size();

    let mut emptyfile = temp_file_manager.create_new_file("empty-", 10, ".bin");
    let char_cursor_positions = crate::bitpacked_kmer_sorting::cursors::init_char_cursor_positions::<B>(&mut emptyfile, sorted_kmers, k, sigma);

    let mut global_cursor = crate::bitpacked_kmer_sorting::cursors::DummyNodeMerger::new(
        &mut emptyfile,
        sorted_kmers,
        k,
    );

    let mut char_cursors = split_global_cursor(&global_cursor, &char_cursor_positions, sigma, k);

    let xs = crate::bitpacked_kmer_sorting::cursors::read_kmers(&mut global_cursor);

    let set_bits = char_cursors.par_iter_mut().enumerate().map(|(c, cursor)| {
        get_set_bits(&xs, cursor, k, c as u8)
    }).flatten().collect::<Vec<usize>>();

    let mut has_predecessor = simple_sds_sbwt::raw_vector::RawVector::new();
    has_predecessor.resize(n, false);

    // TODO there should be a better way to initialize this?
    set_bits.iter().for_each(|idx| {
        has_predecessor.set_bit(*idx, true);
    });

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
