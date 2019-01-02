use fastq_pair::{create_io, delete_empty_fastq, Output, parse_header, parse_read, PartialRead};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Write};
use super::Result;


/// Pair FASTQ files and write out R1/R2 and singleton file
pub fn pair_fastqs(r1_path: &str, r2_path: &str) -> Result<Output> {
    let mut io = create_io(&r1_path, &r2_path)?;
    let mut map = index_read(&mut io.in_read1)?;
    while let Some(read) = parse_read(&mut io.in_read2) {
        // Check if header is in hashmap
        let header = parse_header(&read.header)?;
        if map.contains_key(&header) {
            // Write to BufWriters and remove from HashMap
            let r1 = &map[&header];
            write!(&mut io.out_read1, "{}.1\n{}+\n{}", &header, r1.seq, r1.qscore)?;
            write!(&mut io.out_read2, "{}.2\n{}+\n{}", &header, read.seq, read.qscore)?;
            map.remove(&header);
        } else {
            // Else: Write out R2 to singleton file
            write!(&mut io.out_single, "{}.2\n{}+\n{}", &header, read.seq, read.qscore)?;
        }
    }
    // Write out remainder of singletons left in R1
    for key in map.keys() {
        let r1 = &map[key];
        write!(&mut io.out_single, "{}.1\n{}+\n{}", &key, r1.seq, r1.qscore)?;
        // TODO: Add map.remove() here?
    }
    // Flush output for empty fastq check
    &mut io.out_single.flush()?;

    Ok(Output {
        r1_out_path: io.r1_out_path,
        r2_out_path: io.r2_out_path,
        singleton_path: delete_empty_fastq(&io.singleton_path),
    })
}


/// Create a HashMap associating the unique component of a header to it's
/// sequence and quality score.
pub fn index_read(in_read: &mut BufReader<File>) -> Result<HashMap<String, PartialRead>> {
    let mut map = HashMap::new();
    while let Some(read) = parse_read(in_read) {
        let header = parse_header(&read.header)?;
        map.insert(header, PartialRead { seq: read.seq, qscore: read.qscore });
    }
    Ok(map)
}


#[cfg(test)]
mod tests {
    use std::fs::copy;
    use std::path::Path;
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_store_read1() {
        let file = File::open("data/ncbi_1_paired.fastq").unwrap();
        let mut handle = BufReader::new(file);
        let map = index_read(&mut handle).unwrap();
        let read = &map["@SRR3380692.1"];
        let seq = "ATTGTNTTATTCTATAAAACATTTCAAACCTAGTTAGAGATTTGTAATCAAA\
                    AAACATTTGCGCAGAAAGCAGCACTTAGGGCTGCCTGTTCTATACCCTA\n";
        let qscore = "@@@DD#4AFHHHHJJJJIJJJJJJJJJJJJIIJHGJJIJJJIJJGHGIIJ\
                    JJJIJJJJJJJJIJJHHHFFFFFEEEEEDDDDDDDDDDDCCDEEEFDCDDC\n";
        assert_eq!(seq.to_string(), read.seq);
        assert_eq!(qscore.to_string(), read.qscore);
    }

    #[test]
    fn test_write_pairs() {
        let tmpdir = tempdir().unwrap();
        let tmppath = tmpdir.path();
        let r1_path = tmppath.join("ncbi_1_shuffled.fastq");
        copy("data/ncbi_1_shuffled.fastq", &r1_path).unwrap();
        let input2 = tmppath.join("ncbi_2_shuffled.fastq");
        copy("data/ncbi_2_shuffled.fastq", &input2).unwrap();
        pair_fastqs(r1_path.to_str().unwrap(), input2.to_str().unwrap()).unwrap();
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
    }
}
