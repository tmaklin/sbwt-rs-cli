use std::borrow::BorrowMut;
use std::io::{BufWriter, Write};
use std::io::Seek;
use std::io::SeekFrom;

use super::kmer::LongKmer;
use crate::tempfile::TempFileManager;
use crate::tempfile::TempFile;
use simple_sds_sbwt::raw_vector::*;
use rayon::prelude::*;

use crate::bitpacked_kmer_sorting::cursors::rewind_reader;

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

    global_cursor.dummy_reader.file.rewind();
    global_cursor.nondummy_reader.file.rewind();
    global_cursor.dummy_kmer = DummyNodeMerger::read_from_dummy_reader(global_cursor.dummy_reader);
    global_cursor.nondummy_kmer = DummyNodeMerger::read_from_non_dummy_reader(global_cursor.nondummy_reader, k);
    global_cursor.dummy_position = 0;
    global_cursor.nondummy_position = 0;


    for c in 0..(sigma as u8) {
        let dummy_pos = char_cursors[c as usize].0.1;
        let nondummy_pos = char_cursors[c as usize].1.1;

        global_cursor.dummy_reader.file.set_position(char_cursors[c as usize].0.0);
        global_cursor.nondummy_reader.file.set_position(char_cursors[c as usize].1.0);
        global_cursor.dummy_kmer = DummyNodeMerger::read_from_dummy_reader(global_cursor.dummy_reader);
        global_cursor.nondummy_kmer = DummyNodeMerger::read_from_non_dummy_reader(global_cursor.nondummy_reader, k);
        global_cursor.dummy_position = dummy_pos as usize;
        global_cursor.nondummy_position = nondummy_pos as usize;

        xs.iter().for_each(|(x, _)| {
            let xc = x.set_from_left(k-1, 0).right_shift(1).set_from_left(0, c);

            while global_cursor.peek().is_some(){
                match global_cursor.peek().unwrap().0.cmp(&xc) {
                    std::cmp::Ordering::Greater => break,
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

    // Rewind cursors
    emptyfile.file.seek(SeekFrom::Start(0)).unwrap();
    sorted_kmers.file.seek(SeekFrom::Start(0)).unwrap();

    let mut global_cursor = crate::bitpacked_kmer_sorting::cursors::DummyNodeMerger::new(
        &mut emptyfile,
        sorted_kmers,
        k,
    ); // New global cursor

    // Todo: stream to memory and sort there
    let mut required_dummies = Vec::<(LongKmer::<B>, u8)>::new(); // Pairs (data, length)

    while let Some((x, _)) = global_cursor.peek(){
        if !has_predecessor.bit(global_cursor.nondummy_position()){
            let mut prefix = x;
            for i in 0..k {
                let len = k - i - 1;
                prefix = prefix.left_shift(1);
                required_dummies.push((prefix, len as u8));
            }
        }
        global_cursor.next(); // Advance
    }

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
