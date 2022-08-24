// a set of common utilities for all StanEx subcommands

use path_abs::PathFile;
use std::ffi::OsStr;
use std::process::Command;

// create an absolute file path from a relative file path
// (file must already exist)
pub fn absolute_filepath_checked(relative: &str) -> PathFile {
    match PathFile::new(relative) {
        Ok(absolute) => absolute,
        Err(x) => {
            eprintln!("File does not exist: {}", relative);
            panic!("{}", x);
        }
    }
}

// create an absolute file path from a relative file path
// (file may not already exist)
pub fn absolute_filepath_unchecked(relative: &str) -> PathFile {
    match PathFile::create(relative) {
        Ok(absolute) => absolute,
        Err(x) => {
            eprintln!("An error occurred creating the file: {}", relative);
            panic!("{}", x);
        }
    }
}

// enum struct to represent both single-end and paired-ends reads files
#[derive(Clone)]
pub enum Reads {
    SingleEnd(String),
    PairedEnds(String, String),
}

// creates a bwa index if one does not already exist
// by default, bwa index will create a new file with name equal to the original file name + ".bwt"
pub fn bwa_index_if_required(ref_name: &str) {
    // first, create the absolute filepath from the relative filepath (but throw an error if it doesn't exist)
    let ref_path: PathFile = absolute_filepath_checked(ref_name);
    let ref_path_os_str: &OsStr = ref_path.as_ref();
    let ref_path_str: &str = ref_path_os_str.to_str().unwrap();

    // now check if the BWA index already exists
    let bwa_index_path_str = format!("{}.{}", ref_path_str, "bwt");
    if let Ok(_) = PathFile::new(bwa_index_path_str) {
        println!("BWA index already exists");
        return;
    }

    // Now that we know that we have to index:
    println!("Waiting for bwa index...");
    let mut child_proc = Command::new("bwa")
        .args(&["index", ref_path_str])
        .spawn()
        .unwrap();
    let _result = child_proc.wait().unwrap();
    println!("BWA index complete");
}

// does an alignment using BWA MEM
pub fn bwa_mem_align(ref_name: &str, reads_names: &Reads, result_file: &str, bwa_threads: u16) {
    // first, create the absolute filepaths from the relative filepaths of the ref and reads (throw an error if they don't exist)
    let ref_path: PathFile = absolute_filepath_checked(ref_name);
    let ref_path_os_str: &OsStr = ref_path.as_ref();
    let ref_path_str: &str = ref_path_os_str.to_str().unwrap();

    let absolute_reads: Reads = match reads_names {
        Reads::SingleEnd(filename) => {
            let file_path: PathFile = absolute_filepath_checked(filename);
            let file_path_os_str: &OsStr = file_path.as_ref();
            let file_path_str: &str = file_path_os_str.to_str().unwrap();
            Reads::SingleEnd(file_path_str.to_owned())
        }
        Reads::PairedEnds(file1, file2) => {
            let file1_path: PathFile = absolute_filepath_checked(file1);
            let file1_path_os_str: &OsStr = file1_path.as_ref();
            let file1_path_str: &str = file1_path_os_str.to_str().unwrap();

            let file2_path: PathFile = absolute_filepath_checked(file2);
            let file2_path_os_str: &OsStr = file2_path.as_ref();
            let file2_path_str: &str = file2_path_os_str.to_str().unwrap();

            Reads::PairedEnds(file1_path_str.to_owned(), file2_path_str.to_owned())
        }
    };

    // now do the alignment and store in the result file
    println!("Waiting for bwa mem...");
    match absolute_reads {
        Reads::SingleEnd(filepath) => {
            println!(
                "bwa mem -t {} -o {} {} {}",
                &bwa_threads.to_string()[..],
                &result_file[..],
                ref_path_str,
                &filepath[..]
            );
            let mut child_proc = Command::new("bwa")
                .args(&[
                    "mem",
                    "-t",
                    &bwa_threads.to_string()[..],
                    "-o",
                    &result_file[..],
                    ref_path_str,
                    &filepath[..],
                ])
                .spawn()
                .unwrap();
            let _result = child_proc.wait().unwrap();
        }
        Reads::PairedEnds(file1, file2) => {
            let mut child_proc = Command::new("bwa")
                .args(&[
                    "mem",
                    "-t",
                    &bwa_threads.to_string()[..],
                    "-o",
                    &result_file[..],
                    ref_path_str,
                    &file1[..],
                    &file2[..],
                ])
                .spawn()
                .unwrap();
            let _result = child_proc.wait().unwrap();
        }
    }
    println!("Alignment complete");
}
