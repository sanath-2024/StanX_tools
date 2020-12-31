use path_abs::PathDir;

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::str::FromStr;

use crate::utils;

// tile a set of artificial 150-bp-wide artificial "reads" across the reference genome
// this will then be used by the TE mapper to look for transposons to find
// all possible transposons within the reference genome
// (from the list of possible transposon sequences that we are looking for)

pub fn tile_ref(ref_path: &str, output_dir: &str) {
    let _ref_path_checked = utils::absolute_filepath_checked(ref_path);
    let _output_dir_unchecked = PathDir::create(output_dir);
    let output_path = format!("{}/{}", output_dir, "tiled_ref.fastq");
    let lines = BufReader::new(File::open(ref_path).unwrap())
        .lines()
        .map(|l| l.unwrap());
    let mut writer = BufWriter::new(File::create(output_path).unwrap());
    // store current state information
    // such as chromosome name & length, original position,
    // current read number & name, and char buffer
    let mut chrom = String::new();
    let mut chrom_length: u64 = 0;
    let mut original_pos: u64 = 1;
    let mut read_num: u64 = 1;
    let mut read_name;
    let mut buffer = String::new();
    for line in lines {
        // FASTA header line
        if line.chars().nth(0) == Some('>') {
            let fields: Vec<&str> = (&line[..]).split(" ").collect();
            chrom = fields[0][1..].to_owned();
            chrom_length = FromStr::from_str(
                fields[2]
                    .split("..")
                    .nth(1)
                    .unwrap()
                    .split(";")
                    .nth(0)
                    .unwrap(),
            )
            .unwrap();
            original_pos = 1;
            read_num = 1;
            buffer = String::new();
        } else {
            // process the line char by char
            for nt in line.chars() {
                // print out status every 1,000,000 nts processed
                if original_pos % 1_000_000 == 0 {
                    println!("Processing chromosome {}: position {}", chrom, original_pos);
                }
                // print out status when done processing each chromosome
                if original_pos == chrom_length {
                    println!("Done processing chromosome {}", chrom);
                }
                // add the last character to the buffer
                buffer.push(nt);
                // if the buffer is over-full, remove the first character
                if buffer.len() == 151 {
                    buffer.remove(0);
                }
                // generate fastq read if the buffer is full
                // (use ~ for quality score to indicate highest quality)
                if buffer.len() == 150 {
                    read_name = format!("{}_Read_{}", chrom, read_num);
                    let mut quality_bytes: Vec<u8> = Vec::with_capacity(buffer.len());
                    for _ in 0..buffer.len() {
                        quality_bytes.push(0x7e);
                    }
                    let quality_str = String::from_utf8(quality_bytes).unwrap();
                    let fastq = format!("@{}\n{}\n+\n{}\n", read_name, buffer, quality_str);
                    writer.write(fastq.as_bytes()).unwrap();
                    // update counter and buffer
                    read_num += 1;
                }
                // update counter
                original_pos += 1;
            }
        }
    }
}
