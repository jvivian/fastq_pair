use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufWriter};
use std::io::BufReader;
use failure::{Fallible, ResultExt};
use std::path::Path;
pub type Result<T> = Fallible<T>;

/// A PartialRead which excludes the header to avoid redundancy
/// when storing reads in a HashMap
#[derive(Debug)]
pub struct PartialRead {
    pub seq: String,
    pub qscore: String,
}


// TODO: Make a comprehensive version that also stores descriptor (remainder of header)
/// Structure to hold a single read from a FASTQ file
#[derive(Debug)]
pub struct Read {
    pub header: String,
    pub seq: String,
    pub qscore: String,
}

impl Read {
    fn new() -> Read {
        Read {
            header: String::new(),
            seq: String::new(),
            qscore: String::new(),
        }
    }
}

impl fmt::Display for Read {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n", self.header.trim())?;
        write!(f, "{}\n", self.seq.trim())?;
        write!(f, "+\n")?;
        write!(f, "{}\n", self.qscore.trim())
    }
}


/// Parses a header and returns its unique component
pub fn parse_header(header: &str) -> Result<String> {
    let header_vec: Vec<&str> = header.split_whitespace().collect();
    let uniq = header_vec[0][0..header_vec[0].len() - 2].to_string();
    Ok(uniq)
}


/// Deletes empty FASTQ by parsing read to see if it's reached EOF
pub fn delete_empty_fastq(file_path: &str) -> Result<()> {
    let handle = File::open(file_path)?;
    let mut singleton_reader = BufReader::new(handle);
    match parse_read(&mut singleton_reader) {
        Some(_) => {}
        None => std::fs::remove_file(file_path)?,
    };
    Ok(())
}


/// Parses Read struct from BufReader object
pub fn parse_read(file: &mut impl BufRead) -> Option<Read> {
    let mut sep = String::new();
    let mut read = Read::new();
    // read_line returns a Result<u32> of bytes of the line. EOF is length zero, so we'll break
    if file.read_line(&mut read.header).ok()? == 0 { return None; }
    file.read_line(&mut read.seq).ok()?;
    file.read_line(&mut sep).ok()?;
    file.read_line(&mut read.qscore).ok()?;
    Some(read)
}

/// Contains all Read/Write objects
pub struct IO {
    pub in_read1: BufReader<File>,
    pub in_read2: BufReader<File>,
    pub out_read1: BufWriter<File>,
    pub out_read2: BufWriter<File>,
    pub out_single: BufWriter<File>,
    pub singleton_path: String,
}

/// Create all IO objects for reading and writing
pub fn create_io(r1_path: &str, r2_path: &str) -> Result<IO> {
    let parent = Path::new(r1_path).parent().unwrap();
    let r1_out_path = parent.join("R1_paired.fastq").to_str().unwrap().to_string();
    let r2_out_path = parent.join("R2_paired.fastq").to_str().unwrap().to_string();
    let singleton_path = parent.join("Singletons.fastq").to_str().unwrap().to_string();
    // Readers
    let r1_handle = File::open(&r1_path).context("Can't open read1 file")?;
    let r2_handle = File::open(&r2_path).context("Can't open read2 file")?;
    let r1_reader = BufReader::new(r1_handle);
    let r2_reader = BufReader::new(r2_handle);
    // Writers
    let r1_out_handle = File::create(&r1_out_path).context("Can't create read1 output file")?;
    let r2_out_handle = File::create(&r2_out_path).context("Can't create read2 output file")?;
    let singleton_handle = File::create(&singleton_path).context("Can't create singleton output file")?;
    let r1_writer = BufWriter::new(r1_out_handle);
    let r2_writer = BufWriter::new(r2_out_handle);
    let singleton_writer = BufWriter::new(singleton_handle);
    Ok(IO {
        in_read1: r1_reader,
        in_read2: r2_reader,
        out_read1: r1_writer,
        out_read2: r2_writer,
        out_single: singleton_writer,
        singleton_path: singleton_path.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use super::*;

    #[test]
    fn test_parse_header() {
        let header = "@foo:bar:UUID.1 extra:stuff";
        assert_eq!("@foo:bar:UUID".to_string(), parse_header(header).unwrap());
    }

    #[test]
    fn test_parse_read() {
        let handle = File::open("data/R1.fastq").unwrap();
        let mut reader = BufReader::new(handle);
        let header = "@SRR3380692.1.1 1 length=101\n";
        let seq = "ATTGTNTTATTCTATAAAACATTTCAAACCTAGTTAGAGATTTGTAATCAAA\
                    AAACATTTGCGCAGAAAGCAGCACTTAGGGCTGCCTGTTCTATACCCTA\n";
        let qscore = "@@@DD#4AFHHHHJJJJIJJJJJJJJJJJJIIJHGJJIJJJIJJGHGII\
                    JJJJIJJJJJJJJIJJHHHFFFFFEEEEEDDDDDDDDDDDCCDEEEFDCDDC\n";
        let read = parse_read(&mut reader).unwrap();
        assert_eq!(header.to_string(), read.header);
        assert_eq!(seq.to_string(), read.seq);
        assert_eq!(qscore.to_string(), read.qscore);
    }
}
