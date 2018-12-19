use clap;
use clap::{App, Arg, ArgMatches};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Result};

//use std::process::Command;


// Read Struct
#[derive(Debug, Clone)]
struct PartialRead {
    seq: String,
    qscore: String,
}

impl PartialRead {
    fn new() -> PartialRead {
        PartialRead {
            seq: String::new(),
            qscore: String::new(),
        }
    }
}

// Hold a complete Read
struct Read {
    header: String,
    seq: String,
    qscore: String,
}


fn store_read1(fastq_path: &str) -> Result<HashMap<String, PartialRead>> {
    let mut map = HashMap::new();
    let handle = File::open(fastq_path)?;
    let mut file = BufReader::new(handle);
    loop {
        // Variables for read/hashmap.  I could declare above and run `.clear()` here, but it's more code
        let mut header = String::new();
        let mut sep = String::new();
        let mut read = PartialRead::new();
        // read_line returns a Result<u32> of bytes of the line. EOF is length zero, so we'll break
        if file.read_line(&mut header)? == 0 { break; };
        file.read_line(&mut read.seq)?;
        file.read_line(&mut sep)?;
        file.read_line(&mut read.qscore)?;
        // Parse header for unique component
        let header_vec: Vec<&str> = header.split_whitespace().collect();
        let uniq = header_vec[0].to_string()[0..header_vec[0].len() - 2].to_string(); // Admittedly gross, but I can't use trim_end_matches due to .1.1 repeats
        map.insert(uniq, read.clone());
    }
    Ok(map)
}


fn get_matches() -> ArgMatches<'static> {
    let matches = App::new("fastq_pair")
        .version("1.0")
        .author("John Vivian <jtvivian@gmail.com>")
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
                .required(false)
                .help("Path to Read1 FASTQ")
                .takes_value(true))
        .get_matches();

//    println!("Arguments\n{:#?}", matches);
    matches
}

fn main() {
    // Argument parsing
    let matches = get_matches();

    // Unwrap is safe here due to r1/r2 being required arguments
    let r1_path = matches.value_of("r1").unwrap();
    //let r2_path = matches.value_of("r2").unwrap();

    // Pair fastqs
    let reads = store_read1(&r1_path)?;

    println!("{:?}", reads);
}
