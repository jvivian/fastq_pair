use fastq_pair::{parse_header, parse_read, PartialRead};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Result};
use std::io::Write;
use std::path::Path;

/// Contains all Read/Write objects
struct IO {
    in_read1: BufReader<File>,
    in_read2: BufReader<File>,
    out_read1: BufWriter<File>,
    out_read2: BufWriter<File>,
    out_single: BufWriter<File>,
    singleton_path: String,
}

/// Create all IO objects for reading and writing
fn create_io(r1_path: &str, r2_path: &str) -> Result<IO> {
    let parent = Path::new(r1_path).parent().unwrap();
    let r1_out_path = parent.join("R1_paired.fastq").to_str().unwrap().to_string();
    let r2_out_path = parent.join("R2_paired.fastq").to_str().unwrap().to_string();
    let singleton_path = parent.join("Singletons.fastq").to_str().unwrap().to_string();
    // Readers
    let r1_handle = File::open(&r1_path)?;
    let r2_handle = File::open(&r2_path)?;
    let r1_reader = BufReader::new(r1_handle);
    let r2_reader = BufReader::new(r2_handle);
    // Writers
    let r1_out_handle = File::create(&r1_out_path)?;
    let r2_out_handle = File::create(&r2_out_path)?;
    let singleton_handle = File::create(&singleton_path)?;
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


/// Writes out paired reads two FASTQ files
fn write_read(header: &str,
              w1: &mut BufWriter<File>,
              w2: &mut BufWriter<File>,
              map1: &mut HashMap<String, PartialRead>,
              map2: &mut HashMap<String, PartialRead>) -> Result<()> {
    let r1 = &map1[header];
    let r2 = &map2[header];
    write!(w1, "{}.1\n{}+\n{}", header, r1.seq, r1.qscore)?;
    write!(w2, "{}.2\n{}+\n{}", header, r2.seq, r2.qscore)?;
    &map1.remove(header);
    &map2.remove(header);
    Ok(())
}

fn delete_empty_fastq(file_path: &str) -> Result<()> {
    let handle = File::open(file_path)?;
    let mut singleton_reader = BufReader::new(handle);
    match parse_read(&mut singleton_reader) {
        Some(_) => {}
        None => std::fs::remove_file(file_path)?,
    };
    Ok(())
}

/// Pair two FASTQ files by iterating over both files simultaneously.
/// Should be much more memory efficient than "store_read" method if
/// files are mostly paired
pub fn pair_fastqs(r1_path: &str, r2_path: &str) -> Result<()> {
    let mut io = create_io(r1_path, r2_path)?;
    let mut map1 = HashMap::new();
    let mut map2 = HashMap::new();

    // TODO: Iterate until both reads are empty instead of just R1 -- maybe just use a loop with a break
    while let Some(read1) = parse_read(&mut io.in_read1) {
        let read2 = parse_read(&mut io.in_read2).expect("Failed to read R2");
        let header1 = parse_header(&read1.header)?;
        let header2 = parse_header(&read2.header)?;

        // TODO: Maybe write out read and skip this if header1 == header2 to avoid overhead?
        map1.insert(header1.clone(), PartialRead { seq: read1.seq, qscore: read1.qscore });
        map2.insert(header2.clone(), PartialRead { seq: read2.seq, qscore: read2.qscore });

        // Write reads if header is found in either HashMap
        if map2.contains_key(&header1) {
            write_read(&header1, &mut io.out_read1, &mut io.out_read2, &mut map1, &mut map2)?;
        }
        if map1.contains_key(&header2) {
            write_read(&header2, &mut io.out_read1, &mut io.out_read2, &mut map1, &mut map2)?;
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

    // Delete Singleton file if empty
    io.out_single.flush()?;
    delete_empty_fastq(&io.singleton_path)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_fastqs() {
        pair_fastqs("data/ncbi_1_shuffled.fastq", "data/ncbi_2_shuffled.fastq").unwrap();
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
            assert_eq!(parse_header(&r1.header).unwrap(),
                       parse_header(&r2.header).unwrap());
        }
        // Cleanup
        for output in &outputs {
            std::fs::remove_file(&output).unwrap();
        }
    }
}
