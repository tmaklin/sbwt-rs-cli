use std::io::Cursor;

/// A simple struct to manage temporary files.
/// All files will be created in the directory specified in the struct.
pub struct TempFileManager {
}

/// A struct encapsulating a std::fs::File and its path.
/// When the TempFile is dropped, the file is deleted.
#[derive(Debug, Default)]
pub struct TempFile {
    pub file: Cursor<Vec<u8>>,
}

impl Drop for TempFile {
    fn drop(&mut self) {
	self.file.get_mut().clear();
    }
}

impl std::io::Write for TempFile {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
	self.file.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
	self.file.flush()
    }
}

impl std::io::Read for TempFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
	self.file.read(buf)
    }
}

impl TempFile {
    pub fn avail_in(&self) -> usize {
	return self.file.get_ref().len();
    }
}

impl TempFileManager {

    pub fn new() -> Self {
        Self{}
    }

    // Creates a new temp file with filename with format:
    // self.directory / prefix + random_infix + suffix.
    // &self is taken as mutable because we use an RNG to generate the infix. 
    pub fn create_new_file(&mut self, _prefix: &str, _infix_length: usize, _suffix: &str) -> TempFile {
        TempFile {file: Cursor::new(Vec::new())}
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test_log::test]
//     fn test_tempfile() {
        
//         // Create a random prefix for our temp files. This is to avoid collisions
//         // with files created for previous runs of this test.
//         let seed = rand::thread_rng().gen();
//         let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
//         let mut prefix = rng.gen_range(0, 1000000000).to_string();
//         prefix.push('-');

//         let mut temp_file_manager = TempFileManager::new(Path::new("/tmp"));

//         let mut files = Vec::<TempFile>::new();
//         let mut paths = Vec::<PathBuf>::new();

//         for _ in 0..26 { // Create filename collisions with very high probablity
//             let temp_file = temp_file_manager.create_new_file(&prefix, 1, ".txt");
//             paths.push(temp_file.path.clone());
//             files.push(temp_file);
//         }

//         for path in paths.iter() {
//             assert!(path.exists());
//         }

//         drop(files); // Should delete all our files

//         for path in paths.iter() {
//             assert!(!path.exists());
//         }
//     }
// }
