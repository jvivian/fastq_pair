use clap;
use clap::{App, Arg, ArgMatches};
use fastq_pair::Result;

mod seek;
mod store_read;
mod iter_both;
mod io;

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
                .possible_values(&["store", "seek", "iter"])
                .default_value("store"))
        .arg(
            Arg::with_name("gzip")
                .long("gzip")
                .required(false))
        .get_matches();

    matches
}

fn main() -> Result<()> {
    // Argument parsing
    let matches = cli();
    // Unwrap is safe here due to all arguments being either required or having defaults
    let mut r1_path = matches.value_of("r1").unwrap().to_string();
    let mut r2_path = matches.value_of("r2").unwrap().to_string();
    let method = matches.value_of("method").unwrap();
    let gzip = matches.is_present("gzip");

    // Check input and uncompress if necessary
    r1_path = io::check_input(&r1_path)?;
    r2_path = io::check_input(&r2_path)?;

    // Pair fastqs
    let output;
    match method {
        "store" => {
            output = store_read::pair_fastqs(&r1_path, &r2_path)?;
        }
        "seek" => {
            output = seek::pair_fastqs(&r1_path, &r2_path)?;
        }
        "iter" => {
            output = iter_both::pair_fastqs(&r1_path, &r2_path)?;
        }
        _ => {
            unreachable!();
        }
    }

    // Compress output if selected
    if gzip == true {
        io::gzip(&output.r1_out_path)?;
        io::gzip(&output.r2_out_path)?;
        if let Some(singleton_path) = &output.singleton_path {
            io::gzip(singleton_path)?;
        }
    }

    Ok(())
}
