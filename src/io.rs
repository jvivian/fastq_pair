use std::path::Path;
use subprocess::{Exec, ExitStatus};
use super::Result;

// TODO: Extend to handle BAM
/// Checks input type and processes into fastq if necessary
pub fn convert_to_fastq(input: &str) -> Result<String> {
    let input_path = Path::new(input);
    let ext = input_path.extension()
        .expect("Failed to get extension")
        .to_str().expect("Failed to parse extension to str");
    match ext {
        "gz" => {
            gzip(input)?;
            return Ok(input_path.file_stem()
                .expect("Failed to get file stem")
                .to_str().expect("Failed to parse file stem to str")
                .to_string());
        }
        _ => return Ok(input.to_string())
    };
}

/// Compresses if uncompressed and vice-versa
pub fn gzip(input: &str) -> Result<ExitStatus> {
    let ext = Path::new(input).extension()
        .expect("Failed to get extension")
        .to_str().expect("Failed to parse extension to str");
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
        convert_to_fastq(&path_gz).unwrap();
        assert_eq!(Path::new(&path).exists(), true);
    }

    #[test]
    fn test_gzip() {
        let tmpdir = tempdir().unwrap();
        let tmppath = tmpdir.path();
        let path = tmppath.join("ncbi_1_paired.fastq");
        copy("data/ncbi_1_shuffled.fastq", &path).unwrap();
        let output = tmppath.join("ncbi_1_paired.fastq.gz");
        gzip(path.to_str().unwrap()).unwrap();
        assert_eq!(Path::new(output.to_str().unwrap()).exists(), true);
    }
}
