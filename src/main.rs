use clap;
use clap::{App, Arg, ArgMatches};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Result};
use std::path::Path;


// PartialRead for HashMap
#[derive(Debug)]
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

// Hold a complete Read which
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
fn write_pairs(r2_path: &str, map: HashMap<String, PartialRead>) -> Result<()> {
    // Open R2 and iterate over reads

    // Check if header is in hashmap

    // If in hashmap, write to BufWriter

    // Else: Write out singleton file

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
        let header_vec: Vec<&str> = read.header.split_whitespace().collect();
        let uniq = header_vec[0].to_string()[0..header_vec[0].len() - 2].to_string();
        map.insert(uniq, PartialRead { seq: read.seq, qscore: read.qscore });
    }
    Ok(map)
}

fn parse_read(file: &mut BufReader<File>) -> Option<Read> {
    let mut sep = String::new();
    let mut read = Read::new();
    // read_line returns a Result<u32> of bytes of the line. EOF is length zero, so we'll break
    if file.read_line(&mut read.header).ok()? == 0 { return None; };
    file.read_line(&mut read.seq).ok()?;
    file.read_line(&mut sep).ok()?;
    file.read_line(&mut read.qscore).ok()?;
    // Parse header for unique component
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
                .required(false)
                .help("Path to Read1 FASTQ")
                .takes_value(true))
        .get_matches();

//    println!("Arguments\n{:#?}", matches);
    matches
}

fn main() -> Result<()> {
    // Argument parsing
    let matches = get_matches();

    // Unwrap is safe here due to r1/r2 being required arguments
    let r1_path = matches.value_of("r1").unwrap();
    //let r2_path = matches.value_of("r2").unwrap();

    // Pair fastqs
    let reads = store_read1(&r1_path)?;

    println!("{:?}", reads["@SRR3380692.1"]);

    Ok(())
}
