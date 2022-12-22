use path_abs::PathFile;
use rayon::prelude::*;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};

use super::first_sam_file;

pub fn select_reads(
    te_aligned_path: &PathFile,
    selected_reads_path: &PathFile,
    only_create_transposon_map: bool,
) -> HashMap<String, u64> {
    // select split-reads from TE alignment
    let mut te_aligned_reader =
        BufReader::with_capacity(65_536, File::open(&te_aligned_path).unwrap());
    let selected_reads_writer_arc = Arc::new(Mutex::new(BufWriter::with_capacity(
        65_536,
        File::create(&selected_reads_path).unwrap(),
    )));

    // read the file line by line
    // don't store lines in an intermediate data structure because that wastes memory
    // store the line number in a mutex for later use
    let line_num_arc = Arc::new(Mutex::new(0));

    // first, get rid of comments (comments in the SAM file start with "@SQ")
    // and ignore the last comment line (starts with "@PG")
    // make a clone because transposons will be put into an Arc and cannot be returned
    let transposons = first_sam_file::read_all_tes_into_map(&mut te_aligned_reader);

    if only_create_transposon_map {
        return transposons;
    }

    let transposons_clone = transposons.clone();

    // next, process the normal reads

    // let transposons be borrowed by other threads
    let transposons_arc = Arc::new(transposons);

    te_aligned_reader.lines().par_bridge().for_each(|line| {
        {
            let mut i = line_num_arc.lock().unwrap();
            *i += 1;
            // print status every 1,000,000 lines
            if *i % 1_000_000 == 0 {
                println!("processing line: {}", i);
            }
        }
        let line = line.expect("Something went wrong - unable to read file");
        let alignment = first_sam_file::read_te_alignment(line, &transposons_arc);
        if let Ok(alignment) = alignment {
            selected_reads_writer_arc
                .lock()
                .unwrap()
                .write(format!("{}\n", alignment).as_bytes())
                .unwrap();
        }
    });
    return transposons_clone;
}
