use std::collections::HashMap;
use std::io::{BufRead, Cursor, Seek, SeekFrom, Write};
use std::error::Error;
use super::parse_read;

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
fn pair_fastqs<R, W>(fastq1: &mut R,
                     fastq2: &mut R,
                     paired1: &mut W,
                     paired2: &mut W,
                     unpaired1: &mut W,
                     unpaired2: &mut W) -> Result<(), Box<Error>>
                       where R: Seek + BufRead,
                             W: Write,
{
    let mut index = index_fastq(fastq1);
    while let Some(read2) = parse_read(fastq2) {
        let trimmed = trim_header(&read2.header).unwrap();
        if let Some(pos1) = index.remove(&trimmed) {
            // Pair found -- output them both.
            write!(paired2, "{}", read2)?;
            fastq1.seek(SeekFrom::Start(pos1))?;
            let read1 = parse_read(fastq1).expect("Couldn't read indexed mate");
            write!(paired1, "{}", read1)?;
        } else {
            // No pair detected.
            write!(unpaired2, "{}", read2)?;
        }
    }

    // All the remaining elements of the index are unpaired. Output
    // them into the unpaired file for 1.
    for pos1 in index.drain().map(|(_k, v)| v) {
        fastq1.seek(SeekFrom::Start(pos1))?;
        let read1 = parse_read(fastq1).expect("Couldn't read unpaired mate");
        write!(unpaired1, "{}", read1)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::from_utf8;

    #[test]
    fn test_index_fastq() {
        let fastq = include_str!("../testdata/ncbi_1_shuffled.fastq");
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
        let input1 = include_str!("../testdata/ncbi_1_shuffled.fastq");
        let input2 = include_str!("../testdata/ncbi_2_shuffled.fastq");
        let mut paired1 = Cursor::new(vec![]);
        let mut paired2 = Cursor::new(vec![]);
        let mut unpaired1 = Cursor::new(vec![]);
        let mut unpaired2 = Cursor::new(vec![]);
        pair_fastqs(&mut Cursor::new(input1), &mut Cursor::new(input2),
                    &mut paired1, &mut paired2,
                    &mut unpaired1, &mut unpaired2).expect("Pairing failed");
        assert_eq!(from_utf8(&paired1.into_inner()).unwrap(),
                   include_str!("../testdata/ncbi_1_paired.fastq"));
        assert_eq!(from_utf8(&paired2.into_inner()).unwrap(),
                   include_str!("../testdata/ncbi_2_paired.fastq"));
        assert_eq!(from_utf8(&unpaired1.into_inner()).unwrap(),
                   include_str!("../testdata/ncbi_1_unpaired.fastq"));
        assert_eq!(from_utf8(&unpaired2.into_inner()).unwrap(),
                   include_str!("../testdata/ncbi_2_unpaired.fastq"));
    }
}
