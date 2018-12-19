use std::fs::File;
use std::io::{Result, BufRead, BufReader};
use std::fmt;

// PartialRead for HashMap
#[derive(Debug)]
pub struct PartialRead {
    pub seq: String,
    pub qscore: String,
}

// Complete Read with header
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




pub fn parse_header(header: &str) -> Result<String> {
    let header_vec: Vec<&str> = header.split_whitespace().collect();
    let uniq = header_vec[0][0..header_vec[0].len() - 2].to_string();
    Ok(uniq)
}


/// Parses Read from BufReader object
pub fn parse_read(file: &mut impl BufRead) -> Option<Read> {
    let mut sep = String::new();
    let mut read = Read::new();
    // read_line returns a Result<u32> of bytes of the line. EOF is length zero, so we'll break
    if file.read_line(&mut read.header).ok()? == 0 { return None; };
    file.read_line(&mut read.seq).ok()?;
    file.read_line(&mut sep).ok()?;
    file.read_line(&mut read.qscore).ok()?;
    Some(read)
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
