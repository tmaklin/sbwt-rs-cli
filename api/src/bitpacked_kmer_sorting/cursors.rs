use std::borrow::BorrowMut;
use std::io::Read;
use std::io::Write;

use simple_sds_sbwt::{raw_vector::AccessRaw};
use std::cmp::min;
use super::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;
use crate::tempfile::TempFile;

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;

pub struct DummyNodeMerger<R: std::io::Read, const B: usize> {
    dummy_reader: R, // Stream of k-mer objects
    nondummy_reader: R, // Stream of pairs (kmer, len)

    dummy_kmer: Option<(LongKmer::<B>, u8)>,
    nondummy_kmer: Option<(LongKmer::<B>, u8)>,

    k: usize,

    dummy_position: usize, // Position of the dummy cursor
    nondummy_position: usize, // Position of the nondummy cursor
}

impl <R: std::io::Read, const B: usize> DummyNodeMerger<R, B> {

    pub fn read_from_dummy_reader(dummy_reader: &mut R) -> Option<(LongKmer::<B>, u8)>{
        let kmer = match LongKmer::<B>::load(dummy_reader){
            Ok(kmer_opt) => {
                match kmer_opt{
                    Some(kmer) => kmer,
                    None => return None, // End of stream
                }   
            },
            Err(e) => panic!("IO error while streaming sorted k-mers: {}", e),
        };

        // Read length
        let mut buf = [0_u8; 1];
        dummy_reader.read_exact(&mut buf).unwrap();
        let len = u8::from_le_bytes(buf);

        Some((kmer, len))
    }

    pub fn read_from_non_dummy_reader(nondummy_reader: &mut R, k: usize) -> Option<(LongKmer::<B>, u8)>{
        match LongKmer::<B>::load(nondummy_reader){
            Ok(kmer_opt) => {
                kmer_opt.map(|kmer| (kmer, k as u8))
            },
            Err(e) => panic!("IO error while streaming sorted k-mers: {}", e),
        }
    }

    pub fn new(mut dummy_reader: R, mut nondummy_reader: R, k: usize) -> Self {
        let dummy_kmer = Self::read_from_dummy_reader(&mut dummy_reader);
        let nondummy_kmer = Self::read_from_non_dummy_reader(&mut nondummy_reader, k);

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position: 0,
            nondummy_position: 0,
        }
    }

    pub fn new_with_initial_positions(mut dummy_reader: R, mut nondummy_reader: R, k: usize, dummy_position: usize, nondummy_position: usize) -> Self {
        let dummy_kmer = Self::read_from_dummy_reader(&mut dummy_reader);
        let nondummy_kmer = Self::read_from_non_dummy_reader(&mut nondummy_reader, k);

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position,
            nondummy_position,
        }
    }

    pub fn peek(&self) -> Option<(LongKmer::<B>, u8)>{
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    Some(dummy_kmer)
                } else {
                    Some(nondummy_kmer)
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn dummy_position(&self) -> usize{
        self.dummy_position
    }

    #[allow(dead_code)]
    pub fn nondummy_position(&self)  -> usize{
        self.nondummy_position
    }
}

impl<const B: usize> Iterator for DummyNodeMerger<&mut TempFile, B> {
    type Item = (LongKmer<B>, u8);

    // Produces pairs (kmer, length)
    fn next(&mut self) -> Option<(LongKmer::<B>, u8)> {
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                self.dummy_kmer = Self::read_from_dummy_reader(&mut self.dummy_reader);
                self.dummy_position += 1;
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                self.nondummy_kmer = Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k);
                self.nondummy_position += 1;
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    self.dummy_kmer = Self::read_from_dummy_reader(&mut self.dummy_reader);
                    self.dummy_position += 1;
                    Some(dummy_kmer)
                } else {
                    self.nondummy_kmer = Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k);
                    self.nondummy_position += 1;
                    Some(nondummy_kmer)
                }
            }
        }
    }

}


impl<const B: usize> Iterator for DummyNodeMerger<std::io::Cursor<Vec<u8>>, B> {
    type Item = (LongKmer<B>, u8);

    // Produces pairs (kmer, length)
    fn next(&mut self) -> Option<(LongKmer::<B>, u8)> {
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                self.dummy_kmer = Self::read_from_dummy_reader(&mut self.dummy_reader);
                self.dummy_position += 1;
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                self.nondummy_kmer = Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k);
                self.nondummy_position += 1;
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    self.dummy_kmer = Self::read_from_dummy_reader(&mut self.dummy_reader);
                    self.dummy_position += 1;
                    Some(dummy_kmer)
                } else {
                    self.nondummy_kmer = Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k);
                    self.nondummy_position += 1;
                    Some(nondummy_kmer)
                }
            }
        }
    }

}

pub fn find_in_dummy<const B: usize>(
    dummy_file: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    c: u8,
) -> (u64, u64) {
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

    (start as u64, start as u64)
}

pub fn find_in_nondummy<const B: usize>(
    nondummy_file: &std::io::Cursor<Vec<LongKmer<B>>>,
    c: u8,
) -> (u64, u64) {
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

    (start as u64, start as u64)

}

// We take in Paths instead of a Files because we need multiple readers to the same files 
pub fn init_char_cursor_positions<const B: usize>(
    kmers: &std::io::Cursor<Vec<LongKmer<B>>>,
    dummies: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    _k: usize,
    sigma: usize
) -> Vec<((u64, u64), (u64, u64))>{
    let mut char_cursor_positions = Vec::<((u64, u64), (u64, u64))>::new();
    for c in 0..(sigma as u8){
        log::trace!("Searching character {}", c);

        let (dummy_reader_pos, dummy_pos) = find_in_dummy::<B>(dummies, c);
        let (nondummy_reader_pos, nondummy_pos) = find_in_nondummy::<B>(kmers, c);

        let cursor = ((dummy_reader_pos, dummy_pos), (nondummy_reader_pos, nondummy_pos));
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
    global_cursor: &mut DummyNodeMerger<&mut std::io::Cursor<Vec<u8>>, B>,
    char_cursor_positions: &Vec<((u64, u64), (u64, u64))>,
    sigma: usize,
    k: usize,
) -> Vec<DummyNodeMerger<std::io::Cursor::<Vec<u8>>, B>> {
    let mut char_cursors = (0..(sigma - 1)).collect::<Vec<usize>>().into_iter().map(|c|{
        DummyNodeMerger::new_with_initial_positions(
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.dummy_reader.get_ref()
                    [(char_cursor_positions[c as usize].0.0 as usize)..(char_cursor_positions[c + 1 as usize].0.0 as usize)].to_vec()
            ),
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.nondummy_reader.get_ref()
                    [(char_cursor_positions[c as usize].1.0 as usize)..(char_cursor_positions[c + 1 as usize].1.0 as usize)].to_vec()
            ),
            k,
            char_cursor_positions[c as usize].0.1 as usize,
            char_cursor_positions[c as usize].1.1 as usize,
        )
    }).collect::<Vec<DummyNodeMerger<std::io::Cursor::<Vec<u8>>, B>>>();
    char_cursors.push(
        DummyNodeMerger::new_with_initial_positions(
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.dummy_reader.get_ref()
                    [(char_cursor_positions[sigma - 1 as usize].0.0 as usize)..(global_cursor.dummy_reader.get_ref().len())].to_vec()
            ),
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.nondummy_reader.get_ref()
                    [(char_cursor_positions[sigma - 1 as usize].1.0 as usize)..(global_cursor.nondummy_reader.get_ref().len())].to_vec()
            ),
            k,
            char_cursor_positions[sigma - 1 as usize].0.1 as usize,
            char_cursor_positions[sigma - 1 as usize].1.1 as usize,
        )
    );
    char_cursors
}

pub fn read_kmers<const B: usize>(
    global_cursor: &mut DummyNodeMerger<&mut std::io::Cursor<Vec<u8>>, B>,
    k: usize,
) -> Vec<(LongKmer::<B>, u8)> {
    let n_kmers = global_cursor.nondummy_reader.get_ref().len() / LongKmer::<B>::byte_size();
    let n_dummies = global_cursor.dummy_reader.get_ref().len() / (LongKmer::<B>::byte_size() + 1);
    let n_merged = n_kmers + n_dummies;

    let dummies = global_cursor.dummy_reader.get_mut().par_chunks(LongKmer::<B>::byte_size() + 1).map(|mut bytes| {
        let kmer = LongKmer::<B>::load(&mut bytes).expect("Valid k-mer").unwrap();
        let mut buf = [0_u8; 1];
        bytes.read_exact(&mut buf).unwrap();
        let len = u8::from_le_bytes(buf);
        (kmer, len)
    }).collect::<Vec<(LongKmer::<B>, u8)>>();

    global_cursor.nondummy_reader.set_position(0);

    let mut dummy_idx = 0;
    let mut prev_kmer = (LongKmer::<B>::load(&mut global_cursor.nondummy_reader).expect("Valid k-mer"), k as u8);
    (0..n_merged).map(|_| {
        // Could implement default for LongKmer and see if using mem:;take is faster
        if dummy_idx >= n_dummies {
            let kmer = prev_kmer;
            prev_kmer = (LongKmer::<B>::load(&mut global_cursor.nondummy_reader).expect("Valid k-mer"), k as u8);
            (kmer.0.unwrap(), kmer.1)
        } else if !prev_kmer.0.is_some() {
            dummy_idx += 1;
            dummies[dummy_idx - 1]
        } else if dummies[dummy_idx] < (prev_kmer.0.unwrap(), prev_kmer.1) {
            dummy_idx += 1;
            dummies[dummy_idx - 1]
        } else {
            let kmer = prev_kmer;
            prev_kmer = (LongKmer::<B>::load(&mut global_cursor.nondummy_reader).expect("Valid k-mer"), k as u8);
            (kmer.0.unwrap(), kmer.1)
        }
    }).collect::<Vec<(LongKmer::<B>, u8)>>()
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    kmers_file: &std::io::Cursor<Vec<LongKmer<B>>>,
    required_dummies: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    n: usize,
    k: usize, 
    sigma: usize,
    build_lcs: bool) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>)
{

    let mut serialized_kmers = std::io::Cursor::new(Vec::new());
    kmers_file.get_ref().iter().for_each(|kmer| {
        kmer.serialize(&mut serialized_kmers).expect("Serialized kmer to TempFile");
    });

    let mut serialized_dummies = std::io::Cursor::new(Vec::new());
    required_dummies.get_ref().iter().for_each(|(kmer, len)| {
        kmer.serialize(&mut serialized_dummies).expect("Serialized kmer to TempFile");
        serialized_dummies.write_all(&[*len]).unwrap();
    });

    let mut char_cursor_positions = init_char_cursor_positions::<B>(kmers_file, required_dummies, k, sigma);
    for c in 0..sigma {
        char_cursor_positions[c].0.0 *= LongKmer::<B>::byte_size() as u64 + 1;
        char_cursor_positions[c].1.0 *= LongKmer::<B>::byte_size() as u64;
    }

    let mut global_cursor = DummyNodeMerger::new(
        &mut serialized_dummies,
        &mut serialized_kmers,
        k,
    );

    let mut rawrows = vec![simple_sds_sbwt::raw_vector::RawVector::with_len(n, false); sigma];

    let kmers = read_kmers::<B>(global_cursor.borrow_mut(), k);
    let mut char_cursors = split_global_cursor(&mut global_cursor, &char_cursor_positions, sigma, k);

    char_cursors.iter_mut().zip(rawrows.iter_mut()).enumerate().par_bridge().for_each(|(c, (cursor, rawrows))|{
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

            while cursor.peek().is_some() && cursor.peek().unwrap() < kmer_c {
                cursor.next();
            }

            if cursor.peek().is_some() && cursor.peek().unwrap() == kmer_c {
                rawrows.set_bit(kmer_idx, true);
                cursor.next();
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

// #[cfg(test)]
// mod tests{

//     use std::io::Write;

//     use super::*;

//     #[test]
//     fn test_init_char_cursors(){
//         let nondummies = [
//             LongKmer::<2>::from_ascii(b"ACGT").unwrap(),
//             LongKmer::<2>::from_ascii(b"AGGT").unwrap(),
//             LongKmer::<2>::from_ascii(b"GGAA").unwrap(),
//             LongKmer::<2>::from_ascii(b"GGGT").unwrap()
//         ];
//         let dummies = [
//             (LongKmer::<2>::from_ascii(b"AAAA").unwrap(),0), // This is actually the empty dummy so it's not in the A-block
//             (LongKmer::<2>::from_ascii(b"AAAA").unwrap(),1),
//             (LongKmer::<2>::from_ascii(b"ACAA").unwrap(),2),
//             (LongKmer::<2>::from_ascii(b"ACAA").unwrap(),3),
//             (LongKmer::<2>::from_ascii(b"GGTT").unwrap(),3),
//         ];

//         let mut temp_file_manager = crate::tempfile::TempFileManager::new();

//         let mut nondummy_file = temp_file_manager.create_new_file("test-", 10, ".nondummy");
//         let mut dummy_file = temp_file_manager.create_new_file("test-", 10, ".dummy");

//         for kmer in nondummies.iter(){
//             kmer.serialize(&mut nondummy_file).unwrap();
//         }
//         for (kmer, len) in dummies.iter(){
//             kmer.serialize(&mut dummy_file).unwrap();
//             let len_byte = *len as u8;
//             dummy_file.write_all(&[len_byte]).unwrap();
//         }

//         // Flush
//         dummy_file.flush().unwrap();
//         nondummy_file.flush().unwrap();

//         let char_cursors = init_char_cursors(&mut dummy_file, &mut nondummy_file, 4, 4);

//         assert_eq!(char_cursors[0].peek(), Some((LongKmer::<2>::from_ascii(b"AAAA").unwrap(), 1))); // A
//         assert_eq!(char_cursors[1].peek(), Some((LongKmer::<2>::from_ascii(b"GGAA").unwrap(), 4))); // C
//         assert_eq!(char_cursors[2].peek(), Some((LongKmer::<2>::from_ascii(b"GGAA").unwrap(), 4))); // G
//         assert_eq!(char_cursors[3].peek(), None); // T

//     }
// }
