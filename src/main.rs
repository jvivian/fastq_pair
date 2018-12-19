use clap;
use clap::{App, Arg, ArgMatches};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Result, Write};
use std::path::Path;


// PartialRead for HashMap
#[derive(Debug)]
struct PartialRead {
    seq: String,
    qscore: String,
}

// Complete Read with header
struct Read {
    header: String,
    seq: String,
    qscore: String,
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

// Give ownership of HashMap so I can pop items as pairs are found and it's not used elsewhere
fn write_pairs(r2_path: &str, map: &mut HashMap<String, PartialRead>) -> Result<()> {
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


fn store_read1(r1_path: &str) -> Result<HashMap<String, PartialRead>> {
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

fn parse_header(header: &str) -> Result<String> {
    let header_vec: Vec<&str> = header.split_whitespace().collect();
    let uniq = header_vec[0][0..header_vec[0].len() - 2].to_string();
    Ok(uniq)
}


/// Parses Read from BufReader object
fn parse_read(file: &mut BufReader<File>) -> Option<Read> {
    let mut sep = String::new();
    let mut read = Read::new();
    // read_line returns a Result<u32> of bytes of the line. EOF is length zero, so we'll break
    if file.read_line(&mut read.header).ok()? == 0 { return None; };
    file.read_line(&mut read.seq).ok()?;
    file.read_line(&mut sep).ok()?;
    file.read_line(&mut read.qscore).ok()?;
    Some(read)
}


fn get_matches() -> ArgMatches<'static> {
    let matches = App::new("fastq_pair")
        .version("1.0")
        .author("John Vivian and Joel ArmStrong")
        .about("Pairs two provided fastq files")
        .arg(
            Arg::with_name("r1")
                .short("1")
                .long("read1")
                .value_name("PATH")
                .required(true)
                .help("Path to Read1 FASTQ")
                .takes_value(true))
        .arg(
            Arg::with_name("r2")
                .short("2")
                .long("read2")
                .value_name("PATH")
                .required(true)
                .help("Path to Read1 FASTQ")
                .takes_value(true))
        .get_matches();

    matches
}

fn main() -> Result<()> {
    // Argument parsing
    let matches = get_matches();

    // Unwrap is safe here due to r1/r2 being required arguments
    let r1_path = matches.value_of("r1").unwrap();
    let r2_path = matches.value_of("r2").unwrap();

    // Pair fastqs
    let mut reads = store_read1(&r1_path)?;
    write_pairs(&r2_path, &mut reads)?;

    Ok(())
}
