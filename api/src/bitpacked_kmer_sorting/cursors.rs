use std::borrow::BorrowMut;
use std::io::{Seek, Read};

use simple_sds_sbwt::{ops::Access, raw_vector::AccessRaw};
use std::io::SeekFrom;
use std::cmp::min;
use super::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred_mut;
use crate::tempfile::TempFile;

use rayon::iter::ParallelIterator;
use rayon::iter::ParallelBridge;

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

impl<const B: usize> Iterator for DummyNodeMerger<&mut std::io::Cursor<Vec<u8>>, B> {
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


pub fn rewind_reader<const B: usize>(reader: &mut DummyNodeMerger<&mut TempFile, B>) {
    match reader.dummy_reader.file.rewind() {
        Ok(_) => {
            // ok
        },
        Err(e) => panic!("Error while rewinding dummy reader cursor: {}", e),
    }

    match reader.nondummy_reader.file.rewind() {
        Ok(_) => {
            // ok
        },
        Err(e) => panic!("Error while rewinding nondummy reader cursor: {}", e),
    }

    reader.dummy_kmer = DummyNodeMerger::read_from_dummy_reader(reader.dummy_reader);
    reader.nondummy_kmer = DummyNodeMerger::read_from_non_dummy_reader(reader.nondummy_reader, reader.k);
    reader.dummy_position = 0;
    reader.nondummy_position = 0;
}

pub fn reset_reader_position<const B: usize>(
    reader: &mut DummyNodeMerger<&mut TempFile, B>,
    dummy_reader_pos: u64,
    nondummy_reader_pos: u64,
    dummy_pos: usize,
    nondummy_pos: usize
) {
        reader.dummy_reader.file.set_position(dummy_reader_pos);
        reader.nondummy_reader.file.set_position(nondummy_reader_pos);
        reader.dummy_kmer = DummyNodeMerger::read_from_dummy_reader(reader.dummy_reader);
        reader.nondummy_kmer = DummyNodeMerger::read_from_non_dummy_reader(reader.nondummy_reader, reader.k);
        reader.dummy_position = dummy_pos;
        reader.nondummy_position = nondummy_pos;
}

// We take in Paths instead of a Files because we need multiple readers to the same files 
pub fn init_char_cursor_positions<const B: usize>(dummy_file: &mut TempFile, nondummy_file: &mut TempFile, _k: usize, sigma: usize)
-> Vec<((u64, u64), (u64, u64))>{
    let mut char_cursor_positions = Vec::<((u64, u64), (u64, u64))>::new();
    for c in 0..(sigma as u8){
        log::trace!("Searching character {}", c);

        let (dummy_reader_pos, dummy_pos) =
        { // Seek in dummies

            let dummy_file_len = dummy_file.avail_in() as usize;
            let dummy_record_len = LongKmer::<B>::byte_size() + 1; // Pairs (kmer, len byte)
            assert_eq!(dummy_file_len % dummy_record_len, 0);
    
            let access_fn = |pos| {
                dummy_file.file.seek(SeekFrom::Start(pos as u64 * dummy_record_len as u64)).unwrap();
                let kmer = LongKmer::<B>::load(&mut dummy_file.file).unwrap().unwrap(); // Should never be none because we know the file length

                // Read the length byte
                let mut len_buf = [0_u8; 1];
                dummy_file.file.read_exact(&mut len_buf).unwrap();
                let len = u8::from_le_bytes(len_buf);
                (kmer, len)
            };

            let pred_fn = |kmer: (LongKmer::<B>,u8)| {
                kmer.1 > 0 && kmer.0.get_from_left(0) >= c
            };

            let start = binary_search_leftmost_that_fulfills_pred_mut(
                access_fn, 
                pred_fn, 
                dummy_file_len / dummy_record_len);

            dummy_file.file.seek(SeekFrom::Start(start as u64 * dummy_record_len as u64)).unwrap();
            (dummy_file.file.position(), start as u64)
        };

        let (nondummy_reader_pos, nondummy_pos) =
        { // Seek in nondummies

            let nondummy_file_len = nondummy_file.avail_in() as usize;
            let nondummy_record_len = LongKmer::<B>::byte_size();
            assert_eq!(nondummy_file_len % nondummy_record_len, 0);
    
            let access_fn = |pos| {
		        nondummy_file.file.seek(SeekFrom::Start(pos as u64 * nondummy_record_len as u64)).unwrap();
                LongKmer::<B>::load(&mut nondummy_file.file).unwrap().unwrap() // Should never be None because we know the file length
            };

            let pred_fn = |kmer: LongKmer::<B>| {
                kmer.get_from_left(0) >= c
            };

            let start = binary_search_leftmost_that_fulfills_pred_mut(
                access_fn, 
                pred_fn, 
                nondummy_file_len / nondummy_record_len);

            nondummy_file.file.seek(SeekFrom::Start(start as u64 * nondummy_record_len as u64)).unwrap();
            (nondummy_file.file.position(), start as u64)
        };

        let cursor = ((dummy_reader_pos, dummy_pos), (nondummy_reader_pos, nondummy_pos));
        char_cursor_positions.push(cursor);
    }

    nondummy_file.file.seek(SeekFrom::Start(0)).unwrap();
    char_cursor_positions

}

// Returns the LCS array
pub fn build_lcs_array<const B: usize>(
    kmers: &Vec<(LongKmer::<B>, u8)>,
    k: usize,
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);
    let n_kmers = kmers.len();
    let bitwidth = 64 - (k as u64 - 1).leading_zeros();
    let mut lcs = simple_sds_sbwt::int_vector::IntVector::with_len(n_kmers, bitwidth as usize, 0).unwrap();

    let mut prev_kmer = LongKmer::<B>::from_ascii(b"").unwrap();
    let mut prev_len = 0_usize;
    kmers.iter().enumerate().for_each(|(kmer_idx, (kmer, len))| {
        if kmer_idx > 0 {
            // The longest common suffix is the longest common prefix of reversed k-mers
            let mut lcs_value = LongKmer::<B>::lcp(&prev_kmer, &kmer);
            lcs_value = min(lcs_value, min(prev_len, *len as usize));
            lcs.set(kmer_idx, lcs_value as u64);
        }
        prev_kmer = *kmer;
        prev_len = *len as usize;
    });
    lcs
}

pub fn split_global_cursor<const B: usize>(
    global_cursor: &DummyNodeMerger<&mut TempFile, B>,
    char_cursor_positions: &Vec<((u64, u64), (u64, u64))>,
    sigma: usize,
    k: usize,
) -> Vec<DummyNodeMerger<std::io::Cursor::<Vec<u8>>, B>> {
    let mut char_cursors = (0..(sigma - 1)).collect::<Vec<usize>>().into_iter().map(|c|{
        DummyNodeMerger::new_with_initial_positions(
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.dummy_reader.file.get_ref()
                    [(char_cursor_positions[c as usize].0.0 as usize)..(char_cursor_positions[c + 1 as usize].0.0 as usize)].to_vec()
            ),
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.nondummy_reader.file.get_ref()
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
                global_cursor.dummy_reader.file.get_ref()
                    [(char_cursor_positions[sigma - 1 as usize].0.0 as usize)..(global_cursor.dummy_reader.file.get_ref().len())].to_vec()
            ),
            std::io::Cursor::<Vec<u8>>::new(
                global_cursor.nondummy_reader.file.get_ref()
                    [(char_cursor_positions[sigma - 1 as usize].1.0 as usize)..(global_cursor.nondummy_reader.file.get_ref().len())].to_vec()
            ),
            k,
            char_cursor_positions[sigma - 1 as usize].0.1 as usize,
            char_cursor_positions[sigma - 1 as usize].1.1 as usize,
        )
    );
    char_cursors
}

pub fn read_kmers<const B: usize>(
    global_cursor: &mut DummyNodeMerger<&mut TempFile, B>,
) -> Vec<(LongKmer::<B>, u8)> {
    global_cursor.map(|x| {
        x
    }).collect()
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    mut global_cursor: DummyNodeMerger<&mut TempFile, B>,
    char_cursor_positions: &Vec<((u64, u64), (u64, u64))>,
    n: usize,
    k: usize, 
    sigma: usize,
    build_lcs: bool) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>)
{
    let mut rawrows = vec![simple_sds_sbwt::raw_vector::RawVector::with_len(n, false); sigma];

    let kmers = read_kmers(global_cursor.borrow_mut());

    let mut char_cursors = split_global_cursor(&global_cursor, char_cursor_positions, sigma, k);

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
