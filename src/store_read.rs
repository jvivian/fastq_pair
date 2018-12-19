use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{Result, BufWriter, BufReader, Write};

use fastq_pair::{PartialRead, parse_read, parse_header};

// Give ownership of HashMap so I can pop items as pairs are found and it's not used elsewhere
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
    let mut singleton_writer = BufWriter::new(singleton_handle);
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
            write!(&mut singleton_writer, "{}.2\n{}\n+\n{}", &header, read.seq, read.qscore)?;
        }
    }

    // Write out remainder of singletons left in R1
    for key in map.keys() {
        let r1 = &map[key];
        write!(&mut r1_writer, "{}.1\n{}+\n{}", &key, r1.seq, r1.qscore)?;
    }

    // If singleton's file is empty, delete
    let s_handle = File::open(&singleton_path)?;
    let mut s_reader = BufReader::new(s_handle);
    match parse_read(&mut s_reader) {
        Some(_) => {}
        None => std::fs::remove_file(singleton_path)?
    };

    Ok(())
}


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