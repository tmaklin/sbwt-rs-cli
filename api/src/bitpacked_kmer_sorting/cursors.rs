use simple_sds_sbwt::{raw_vector::AccessRaw};
use std::cmp::min;
use super::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;

pub fn find_in_dummy<const B: usize>(
    dummy_file: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    c: u8,
) -> u64 {
    let dummy_file_len = dummy_file.get_ref().len();

    let access_fn = |pos| {
        let record = dummy_file.get_ref()[pos];
        record
    };

    let pred_fn = |kmer: (LongKmer::<B>,u8)| {
        kmer.1 > 0 && kmer.0.get_from_left(0) >= c
    };

    let start = binary_search_leftmost_that_fulfills_pred(
        access_fn,
        pred_fn,
        dummy_file_len);

    start as u64
}

pub fn find_in_nondummy<const B: usize>(
    nondummy_file: &std::io::Cursor<Vec<LongKmer<B>>>,
    c: u8,
) -> u64 {
    let nondummy_file_len = nondummy_file.get_ref().len() as usize;

    let access_fn = |pos| {
        nondummy_file.get_ref()[pos as usize]
    };

    let pred_fn = |kmer: LongKmer::<B>| {
        kmer.get_from_left(0) >= c
    };

    let start = binary_search_leftmost_that_fulfills_pred(
        access_fn,
        pred_fn,
        nondummy_file_len);

    start as u64

}

// We take in Paths instead of a Files because we need multiple readers to the same files 
pub fn init_char_cursor_positions<const B: usize>(
    kmers: &std::io::Cursor<Vec<LongKmer<B>>>,
    dummies: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    sigma: usize
) -> Vec<(u64, u64)>{
    let mut char_cursor_positions = Vec::<(u64, u64)>::new();
    for c in 0..(sigma as u8){
        log::trace!("Searching character {}", c);

        let dummy_reader_pos = find_in_dummy::<B>(dummies, c);
        let nondummy_reader_pos = find_in_nondummy::<B>(kmers, c);

        let cursor = (dummy_reader_pos, nondummy_reader_pos);
        char_cursor_positions.push(cursor);
    }

    char_cursor_positions
}

// Returns the LCS array
pub fn build_lcs_array<const B: usize>(
    kmers: &Vec<(LongKmer::<B>, u8)>,
    k: usize,
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);

    std::iter::once(0_usize).chain(kmers.par_windows(2).map(|window| {
        // The longest common suffix is the longest common prefix of reversed k-mers
        let prev_kmer = &window[0].0;
        let prev_len = &window[0].1;
        let kmer = &window[1].0;
        let len = &window[1].1;
        let mut lcs_value = LongKmer::<B>::lcp(&prev_kmer, &kmer);
        lcs_value = min(lcs_value, min(*prev_len as usize, *len as usize));
        lcs_value
    }).collect::<Vec<usize>>().into_iter()).collect::<simple_sds_sbwt::int_vector::IntVector>()
}

pub fn split_global_cursor<const B: usize>(
    kmers: &mut std::io::Cursor<Vec<LongKmer<B>>>,
    dummies: &mut std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    char_cursor_positions: &Vec<(u64, u64)>,
    sigma: usize,
) -> Vec<(std::io::Cursor<Vec<(LongKmer<B>, u8)>>, std::io::Cursor<Vec<LongKmer<B>>>)> {

    let mut char_cursors = (0..(sigma - 1)).map(|c|{
        (std::io::Cursor::new(Vec::from(
            dummies.get_ref()
                [(char_cursor_positions[c as usize].0 as usize)..(char_cursor_positions[c + 1 as usize].0 as usize)].to_vec()
        )),
        std::io::Cursor::new(Vec::from(
            kmers.get_ref()
                [(char_cursor_positions[c as usize].1 as usize)..(char_cursor_positions[c + 1 as usize].1 as usize)].to_vec()
        )))
    }).collect::<Vec<(std::io::Cursor<Vec<(LongKmer<B>, u8)>>, std::io::Cursor<Vec<LongKmer<B>>>)>>();
    char_cursors.push(
        (std::io::Cursor::new(Vec::from(
            dummies.get_ref()
                [(char_cursor_positions[sigma - 1 as usize].0 as usize)..(dummies.get_ref().len())].to_vec()
        )),
        std::io::Cursor::new(Vec::from(
            kmers.get_ref()
                [(char_cursor_positions[sigma - 1 as usize].1 as usize)..(kmers.get_ref().len())].to_vec()
        ))));

    char_cursors
}

pub fn read_kmer_or_dummy<const B: usize>(
    kmers: &mut std::io::Cursor<Vec<LongKmer<B>>>,
    dummies: &mut std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    k: usize,
) -> Option<(LongKmer<B>, u8)> {
    let n_kmers = kmers.get_ref().len();
    let n_dummies = dummies.get_ref().len();
    let kmers_pos = kmers.position() as usize;
    let dummies_pos = dummies.position() as usize;

    if kmers_pos >= n_kmers && dummies_pos >= n_dummies {
        None
    } else if kmers_pos == n_kmers {
        dummies.set_position(dummies_pos as u64 + 1);
        Some(dummies.get_ref()[dummies_pos])
    } else if dummies_pos == n_dummies {
        kmers.set_position(kmers_pos as u64 + 1);
        Some((kmers.get_ref()[kmers_pos], k as u8))
    } else {
        if dummies.get_ref()[dummies_pos] < (kmers.get_ref()[kmers_pos], k as u8) {
            dummies.set_position(dummies_pos as u64 + 1);
            Some(dummies.get_ref()[dummies_pos])
        } else {
            kmers.set_position(kmers_pos as u64 + 1);
            Some((kmers.get_ref()[kmers_pos], k as u8))
        }
    }
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    kmers_file: &mut std::io::Cursor<Vec<LongKmer<B>>>,
    required_dummies: &mut std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    n: usize,
    k: usize, 
    sigma: usize,
    build_lcs: bool) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>
) {

    let mut kmers: Vec<(LongKmer::<B>, u8)> = Vec::with_capacity(kmers_file.get_ref().len() + required_dummies.get_ref().len());
    while let Some(kmer) = read_kmer_or_dummy(kmers_file, required_dummies, k) {
        kmers.push(kmer);
    }

    let char_cursor_positions = init_char_cursor_positions::<B>(kmers_file, required_dummies, sigma);

    let mut rawrows = vec![simple_sds_sbwt::raw_vector::RawVector::with_len(n, false); sigma];
    let mut char_cursors = split_global_cursor(kmers_file, required_dummies, &char_cursor_positions, sigma);
    char_cursors.iter_mut().zip(rawrows.iter_mut()).enumerate().par_bridge().for_each(|(c, (cursor, rawrows))|{
        let mut pointed_kmer = read_kmer_or_dummy(&mut cursor.1, &mut cursor.0, k);
        kmers.iter().enumerate().for_each(|(kmer_idx, (kmer, len))| {
            let kmer_c = if *len as usize == k {
                (
                    kmer
                        .set_from_left(k - 1, 0)
                        .right_shift(1)
                        .set_from_left(0, c as u8),
                    k as u8,
                )
            } else {
                (kmer.right_shift(1).set_from_left(0, c as u8), len + 1) // Dummy
            };

            while pointed_kmer.is_some() && pointed_kmer.unwrap() < kmer_c {
                pointed_kmer = read_kmer_or_dummy(&mut cursor.1, &mut cursor.0, k);
            }

            if pointed_kmer.is_some() && pointed_kmer.unwrap() == kmer_c {
                rawrows.set_bit(kmer_idx, true);
                pointed_kmer = read_kmer_or_dummy(&mut cursor.1, &mut cursor.0, k);
            }
        });
    });

    let lcs = if build_lcs {
        // LCS values are between 0 and k-1
        Some(build_lcs_array(&kmers, k))
    } else {
        None
    };

    (rawrows, lcs)

}
