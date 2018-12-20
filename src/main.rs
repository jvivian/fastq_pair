use clap;
use clap::{App, Arg, ArgMatches};
use std::io::Result;

mod seek;
mod store_read;

fn cli() -> ArgMatches<'static> {
    let matches = App::new("fastq_pair")
        .version("1.0")
        .author("John Vivian and Joel Armstrong")
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
        .arg(
            Arg::with_name("method")
                .required(false)
                .takes_value(true)
                .possible_values(&["store", "seek"])
                .default_value("store"))
        .get_matches();

    matches
}

fn main() -> Result<()> {
    // Argument parsing
    let matches = cli();
    // Unwrap is safe here due to all arguments being either required or having defaults
    let r1_path = matches.value_of("r1").unwrap();
    let r2_path = matches.value_of("r2").unwrap();
    let method = matches.value_of("method").unwrap();

    // Pair fastqs
    match method {
        "store" => {
            let mut reads = store_read::store_read1(&r1_path)?;
            store_read::write_pairs(&r2_path, &mut reads)?;
        }
        "seek" => {
            // FIXME: main returns an io::Error so this should as
            // well, for consistency
            seek::pair_files(&r1_path, &r2_path).expect("Failed to pair");
        }
        _ => {
            unreachable!();
        }
    }
    Ok(())
}
