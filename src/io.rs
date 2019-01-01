use std::path::Path;
use subprocess::{Exec, ExitStatus};
use super::Result;

/// Checks input type and performs appropriate action
/// Can later be extended to handle BAM or other input
pub fn check_input(input: &str) -> Result<String> {
    let input_path = Path::new(input);
    let ext = input_path.extension().unwrap().to_str().unwrap();
    match ext {
        "gz" => {
            gzip(input)?;
            let a = input_path.file_stem().unwrap().to_str().unwrap().to_string();
            return Ok(a);
        }
        _ => return Ok(input.to_string())
    };
}

/// Compresses if uncompressed and vice-versa
pub fn gzip(input: &str) -> Result<ExitStatus> {
    let ext = Path::new(input).extension().unwrap().to_str().unwrap();
    let command;
    if ext == "gz" { command = "gunzip".to_string(); } else { command = "gzip".to_string(); }
    let exit_status = Exec::cmd(command)
        .arg(input)
        .join()?;
    Ok(exit_status)
}


#[cfg(test)]
mod tests {
    use std::fs::copy;
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_check_input() {
        // Prepare tests
        let tmpdir = tempdir().unwrap();
        let tmppath = tmpdir.path();
        let path = tmppath.join("ncbi_1_paired.fastq");
        copy("data/ncbi_1_paired.fastq", &path).unwrap();
        let path_gz = format!("{}.gz", path.to_str().unwrap());
        gzip(path.to_str().unwrap()).unwrap();
        check_input(&path_gz).unwrap();
        assert_eq!(Path::new(&path).exists(), true);
    }

    #[test]
    fn test_gzip() {
        // Prepare tmppaths
        let tmpdir = tempdir().unwrap();
        let tmppath = tmpdir.path();
        let path = tmppath.join("ncbi_1_paired.fastq");
        copy("data/ncbi_1_shuffled.fastq", &path).unwrap();
        //
        let output = tmppath.join("ncbi_1_paired.fastq.gz");
        gzip(path.to_str().unwrap()).unwrap();
        assert_eq!(Path::new(output.to_str().unwrap()).exists(), true);
    }
}
