use path_abs::PathFile;
use std::ffi::OsStr;
use std::process::Command;

use crate::utils;

pub fn download(url: &str, output_file: &str) {
    // first, create the absolute filepath from the relative filepath (create it if it doesn't exist)
    let output_path: PathFile = utils::absolute_filepath_unchecked(output_file);
    let output_path_os_str: &OsStr = output_path.as_ref();
    let output_path_str: &str = output_path_os_str.to_str().unwrap();
    // curl:
    // -L argument is the location
    // -o argument is the output file
    println!("Waiting for cURL command to download file...");
    let mut child_proc = Command::new("curl")
        .args(&["-L", url, "-o", output_path_str])
        .spawn()
        .unwrap();
    let _result = child_proc.wait().unwrap();
    println!("Reference sequence downloaded");
}
