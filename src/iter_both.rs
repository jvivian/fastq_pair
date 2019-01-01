use fastq_pair::{create_io, delete_empty_fastq, Output, parse_header, parse_read, PartialRead};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use super::Result;

/// Writes out paired reads two FASTQ files
fn write_read(header: &str,
              w1: &mut BufWriter<File>,
              w2: &mut BufWriter<File>,
              map1: &mut HashMap<String, PartialRead>,
              map2: &mut HashMap<String, PartialRead>) -> Result<()> {
    let r1 = &map1.remove(header).expect("Failed to remove header from r1 hashamp");
    let r2 = &map2.remove(header).expect("Failed to remove header from r2 hashmap");
    write!(w1, "{}.1\n{}+\n{}", header, r1.seq, r1.qscore)?;
    write!(w2, "{}.2\n{}+\n{}", header, r2.seq, r2.qscore)?;
    Ok(())
}

/// Pair two FASTQ files by iterating over both files simultaneously.
/// Should be much more memory efficient than "store_read" method if
/// files are mostly paired
pub fn pair_fastqs(r1_path: &str, r2_path: &str) -> Result<Output> {
    let mut io = create_io(r1_path, r2_path)?;
    let mut map1 = HashMap::new();
    let mut map2 = HashMap::new();
    let (mut read1_finished, mut read2_finished) = (false, false);
    while !(read1_finished && read2_finished) {
        if let Some(read1) = parse_read(&mut io.in_read1) {
            let header1 = parse_header(&read1.header)?;
            map1.insert(header1.clone(), PartialRead { seq: read1.seq, qscore: read1.qscore });
            if map2.contains_key(&header1) {
                write_read(&header1, &mut io.out_read1, &mut io.out_read2, &mut map1, &mut map2)?;
            }
        } else { read1_finished = true }
        if let Some(read2) = parse_read(&mut io.in_read2) {
            let header2 = parse_header(&read2.header)?;
            map2.insert(header2.clone(), PartialRead { seq: read2.seq, qscore: read2.qscore });
            if map1.contains_key(&header2) {
                write_read(&header2, &mut io.out_read1, &mut io.out_read2, &mut map1, &mut map2)?;
            } else { read2_finished = true }
        }
    }
    // Write out singletons
    for key in map1.keys() {
        let r1 = &map1[key];
        write!(&mut io.out_single, "{}.1\n{}+\n{}", &key, r1.seq, r1.qscore)?;
    }
    for key in map2.keys() {
        let r2 = &map2[key];
        write!(&mut io.out_single, "{}.2\n{}+\n{}", &key, r2.seq, r2.qscore)?;
    }
    // Flush output for empty fastq check
    io.out_single.flush()?;

    Ok(Output {
        r1_out_path: io.r1_out_path,
        r2_out_path: io.r2_out_path,
        singleton_path: delete_empty_fastq(&io.singleton_path),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use std::io::BufReader;
    use tempfile::tempdir;
    use std::fs::copy;

    #[test]
    fn test_pair_fastqs() {
        let tmpdir = tempdir().unwrap();
        let tmppath = tmpdir.path();
        let input1 = tmppath.join("ncbi_1_shuffled.fastq");
        let input2 = tmppath.join("ncbi_2_shuffled.fastq");
        copy("data/ncbi_1_shuffled.fastq", &input1).unwrap();
        copy("data/ncbi_2_shuffled.fastq", &input2).unwrap();
        pair_fastqs(&input1.to_str().unwrap(), &input2.to_str().unwrap()).unwrap();
        // Output exists
        let outputs = [tmppath.join("R1_paired.fastq"),
                       tmppath.join("R2_paired.fastq"),
                       tmppath.join("Singletons.fastq")];
        for output in &outputs {
            assert_eq!(Path::exists(Path::new(output)), true);
        }
        // Assert each header in each file matches
        let h1 = File::open(&outputs[0]).unwrap();
        let h2 = File::open(&outputs[1]).unwrap();
        let mut reader1 = BufReader::new(&h1);
        let mut reader2 = BufReader::new(&h2);
        for _ in 0..4 {
            let r1 = parse_read(&mut reader1).unwrap();
            let r2 = parse_read(&mut reader2).unwrap();
            assert_eq!(parse_header(&r1.header).unwrap(),
                       parse_header(&r2.header).unwrap());
        }
        // Cleanup
        for output in &outputs {
            std::fs::remove_file(&output).unwrap();
        }
    }
}
