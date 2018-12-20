use fastq_pair::{parse_header, parse_read, PartialRead};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Result, Write};
use std::path::Path;


/// Pair FASTQ files and write out R1/R2 and singleton file
pub fn write_pairs(r2_path: &str, map: &mut HashMap<String, PartialRead>) -> Result<()> {
    // Create output paths
    let parent = Path::new(r2_path).parent().unwrap();
    let r1_out_path = parent.join("R1_paired.fastq").to_str().unwrap().to_string(); // This is sort of ridiculous
    let r2_out_path = parent.join("R2_paired.fastq").to_str().unwrap().to_string();
    let singleton_path = parent.join("Singletons.fastq").to_str().unwrap().to_string();
    // Create BufWriter objects for R1/R2
    let r1_out_handle = File::create(r1_out_path)?;
    let r2_out_handle = File::create(r2_out_path)?;
    let singleton_handle = File::create(&singleton_path)?;
    let mut r1_writer = BufWriter::new(r1_out_handle);
    let mut r2_writer = BufWriter::new(r2_out_handle);
    let mut singleton_writer = BufWriter::new(&singleton_handle);
    // Open R2 and iterate over reads
    let r2_handle = File::open(r2_path)?;
    let mut reader = BufReader::new(r2_handle);
    loop {
        let read = match parse_read(&mut reader) {
            Some(i) => i,
            None => break
        };
        // Check if header is in hashmap
        let header = parse_header(&read.header)?;
        if map.contains_key(&header) {
            // Write to BufWriters and remove from HashMap
            let r1 = &map[&header];
            write!(&mut r1_writer, "{}.1\n{}+\n{}", &header, r1.seq, r1.qscore)?;
            write!(&mut r2_writer, "{}.2\n{}+\n{}", &header, read.seq, read.qscore)?;
            map.remove(&header);
        } else {
            // Else: Write out R2 to singleton file
            write!(&mut singleton_writer, "{}.2\n{}+\n{}", &header, read.seq, read.qscore)?;
        }
    }
    // Write out remainder of singletons left in R1
    for key in map.keys() {
        let r1 = &map[key];
        write!(&mut singleton_writer, "{}.1\n{}+\n{}", &key, r1.seq, r1.qscore)?;
    }
    // If singleton's file is empty, delete
    singleton_writer.flush()?;
    let s_handle = File::open(&singleton_path)?;
    let mut singleton_reader = BufReader::new(s_handle);
    match parse_read(&mut singleton_reader) {
        Some(_) => {}
        None => {
            std::fs::remove_file(singleton_path)?
        }
    };
    Ok(())
}


/// Create a HashMap associating the unique component of a header to it's
/// sequence and quality score.
pub fn store_read1(r1_path: &str) -> Result<HashMap<String, PartialRead>> {
    let mut map = HashMap::new();
    let handle = File::open(r1_path)?;
    let mut file = BufReader::new(handle);
    loop {
        let read = match parse_read(&mut file) {
            Some(i) => i,
            None => break
        };
        // Parse header into unique String and store along with PartialRead in HashMap
        let header = parse_header(&read.header)?;
        map.insert(header, PartialRead { seq: read.seq, qscore: read.qscore });
    }
    Ok(map)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_store_read1() {
        let map = store_read1("data/R1.fastq").unwrap();
        let read = &map["@SRR3380692.1"];
        let seq = "ATTGTNTTATTCTATAAAACATTTCAAACCTAGTTAGAGATTTGTAATCAAA\
                    AAACATTTGCGCAGAAAGCAGCACTTAGGGCTGCCTGTTCTATACCCTA\n";
        let qscore = "@@@DD#4AFHHHHJJJJIJJJJJJJJJJJJIIJHGJJIJJJIJJGHGII\
                    JJJJIJJJJJJJJIJJHHHFFFFFEEEEEDDDDDDDDDDDCCDEEEFDCDDC\n";
        assert_eq!(seq.to_string(), read.seq);
        assert_eq!(qscore.to_string(), read.qscore);
    }

    #[test]
    fn test_write_pairs() {
        let mut map = store_read1("data/ncbi_1_shuffled.fastq").unwrap();
        write_pairs("data/ncbi_2_shuffled.fastq", &mut map).unwrap();
        // Output exists
        let outputs = ["data/R1_paired.fastq", "data/R2_paired.fastq", "data/Singletons.fastq"];
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
            assert_eq!(parse_header(&r1.header).unwrap(), parse_header(&r2.header).unwrap());
        }
        // Cleanup
        for output in &outputs {
            std::fs::remove_file(&output).unwrap();
        }
    }
}
