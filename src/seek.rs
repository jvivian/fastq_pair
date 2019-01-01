use fastq_pair::{create_io, delete_empty_fastq, Output, parse_read};
use std::collections::HashMap;
use std::io::{BufRead, Seek, SeekFrom, Write};
use super::Result;

// FIXME: unify separate implementations.
fn get_next_header(input: &mut impl BufRead) -> Option<String> {
    let mut full_header = String::new();
    if input.read_line(&mut full_header).ok()? == 0 {
        // At EOF
        return None;
    }
    let header = trim_header(&full_header);
    // Skip other 4 lines
    let mut _ignored = String::new();
    input.read_line(&mut _ignored).ok()?;
    input.read_line(&mut _ignored).ok()?;
    input.read_line(&mut _ignored).ok()?;
    header
}

/// From a full header like "@SRR3380692.5.2 3 length=101", get the
/// part that should be identical among mates ("SRR3380692.5").
fn trim_header(full_header: &str) -> Option<String> {
    full_header
        .split(' ')
        .next()
        .map(|h| h.chars().skip(1).take(h.len() - 3).collect())
}

/// Convenience function giving the current offset within a file.
fn tell(f: &mut impl Seek) -> u64 {
    f.seek(SeekFrom::Current(0)).expect("Couldn't seek within file")
}

/// Create an index mapping mates to their location within a file.
fn index_fastq<T>(input: &mut T) -> HashMap<String, u64> where T: Seek + BufRead {
    let mut cur_pos = tell(input);
    let mut index = HashMap::new();
    while let Some(header) = get_next_header(input) {
        index.insert(header, cur_pos);
        cur_pos = tell(input);
    }
    index
}

/// Pair input FASTQ files in a low-memory fashion, writing mates to
/// paired1 and paired2 in the same order. Unpaired reads are output
/// to unpaired1 and unpaired2.
pub fn pair_fastqs(path1: &str, path2: &str) -> Result<Output> {
    let mut io = create_io(path1, path2)?;
    let mut index = index_fastq(&mut io.in_read1);
    while let Some(read2) = parse_read(&mut io.in_read2) {
        let trimmed = trim_header(&read2.header).unwrap();
        if let Some(pos1) = index.remove(&trimmed) {
            // Pair found -- output them both.
            write!(&mut io.out_read2, "{}", read2)?;
            &mut io.in_read1.seek(SeekFrom::Start(pos1))?;
            let read1 = parse_read(&mut io.in_read1).expect("Couldn't read indexed mate");
            write!(&mut io.out_read1, "{}", read1)?;
        } else {
            // No pair detected.
            write!(&mut io.out_single, "{}", read2)?;
        }
    }

    // All the remaining elements of the index are unpaired. Output
    // them into the unpaired file for 1.
    for pos1 in index.drain().map(|(_k, v)| v) {
        &mut io.in_read1.seek(SeekFrom::Start(pos1))?;
        let read1 = parse_read(&mut io.in_read1).expect("Couldn't read unpaired mate");
        write!(&mut io.out_single, "{}", read1)?;
    }
    // Flush output for empty fastq check
    &mut io.out_single.flush()?;

    Ok(Output {
        r1_out_path: io.r1_out_path,
        r2_out_path: io.r2_out_path,
        singleton_path: delete_empty_fastq(&io.singleton_path),
    })
}

#[cfg(test)]
mod tests {
    use std::fs::copy;
    use std::fs::File;
    use std::io::Cursor;
    use std::io::Read;
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_index_fastq() {
        let fastq = include_str!("../data/ncbi_1_shuffled.fastq");
        let read_pos = index_fastq(&mut Cursor::new(fastq.as_bytes()));
        assert_eq!(read_pos, vec![
            ("SRR3380692.3".to_string(), 0),
            ("SRR3380692.2".to_string(), 262),
            ("SRR3380692.1".to_string(), 524),
            ("SRR3380692.4".to_string(), 786),
            ("SRR3380692.9".to_string(), 1048),
        ].into_iter().collect());
    }

    #[test]
    fn test_pair_fastqs() {
        let tmpdir = tempdir().unwrap();
        let tmppath = tmpdir.path();

        let r1_path = tmppath.join("ncbi_1_shuffled.fastq");
        let r2_path = tmppath.join("ncbi_2_shuffled.fastq");
        copy("data/ncbi_1_shuffled.fastq", &r1_path).unwrap();
        copy("data/ncbi_2_shuffled.fastq", &r2_path).unwrap();
        let output = pair_fastqs(r1_path.to_str().unwrap(), r2_path.to_str().unwrap()).expect("Pairing failed");

        let (mut paired1, mut paired2, mut unpaired) = (String::new(), String::new(), String::new());
        File::open(&output.r1_out_path).unwrap().read_to_string(&mut paired1).unwrap();
        File::open(&output.r2_out_path).unwrap().read_to_string(&mut paired2).unwrap();
        let singleton_path = output.singleton_path.unwrap();
        File::open(&singleton_path).unwrap().read_to_string(&mut unpaired).unwrap();
        assert_eq!(&paired1, include_str!("../data/ncbi_1_paired.fastq"));
        assert_eq!(&paired2, include_str!("../data/ncbi_2_paired.fastq"));
        assert_eq!(&unpaired, include_str!("../data/ncbi_unpaired.fastq"));
    }
}
