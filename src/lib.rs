use std::fmt;
use std::fs::File;
use std::io::{BufRead, Result};
use std::io::BufReader;

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
