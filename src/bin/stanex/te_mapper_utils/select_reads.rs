use path_abs::PathFile;
use threadpool::ThreadPool;

use std::collections::HashMap;
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};

use super::preprocess_select_reads::{create_transposon, TeAlignment};

pub fn select_reads(
    te_aligned_path: &PathFile,
    selected_reads_path: &PathFile,
    num_threads: i32,
) -> HashMap<String, u64> {
    // select split-reads from TE alignment
    let te_aligned_reader_arc = Arc::new(Mutex::new(BufReader::with_capacity(
        65_536,
        File::open(&te_aligned_path).unwrap(),
    )));
    let mut te_aligned_reader_main = te_aligned_reader_arc.lock().unwrap();
    let selected_reads_writer_arc = Arc::new(Mutex::new(BufWriter::with_capacity(
        65_536,
        File::create(&selected_reads_path).unwrap(),
    )));

    // read the file line by line
    // don't store lines in an intermediate data structure because that wastes memory
    // store the line number in a mutex for later use
    let line_num_arc = Arc::new(Mutex::new(0));
    let line_num_arc_clone_main = Arc::clone(&line_num_arc);
    let mut line_num_main = line_num_arc_clone_main.lock().unwrap();
    let mut te_aligned_read;

    // first, get rid of comments (comments in the SAM file start with "@SQ")
    // and ignore the last comment line (starts with "@PG")
    let mut transposons: HashMap<String, u64> = HashMap::new();
    // need to make a clone in order to return it at the end of the function
    let mut transposons_clone: HashMap<String, u64> = HashMap::new();
    loop {
        te_aligned_read = String::new();
        te_aligned_reader_main
            .read_line(&mut te_aligned_read)
            .unwrap();
        *line_num_main += 1;
        if te_aligned_read.chars().nth(1).unwrap() == 'P' {
            break;
        } else {
            create_transposon(&te_aligned_read, &mut transposons);
            create_transposon(&te_aligned_read, &mut transposons_clone);
        }
    }

    // next, process the normal reads

    // unlock the mutexes
    std::mem::drop(te_aligned_reader_main);
    std::mem::drop(line_num_main);

    // let transposons be borrowed
    let transposons_arc = Arc::new(transposons);

    // create a threadpool with num_threads workers
    let num_workers: usize = if num_threads <= 0 {
        8
    } else {
        num_threads.try_into().unwrap()
    };
    let pool = ThreadPool::with_name("stanex_tools worker".into(), num_workers);
    for _ in 0..num_workers {
        let te_aligned_reader_arc_clone = Arc::clone(&te_aligned_reader_arc);
        let line_num_arc_clone = Arc::clone(&line_num_arc);
        let selected_reads_writer_arc_clone = Arc::clone(&selected_reads_writer_arc);
        let transposons_arc_clone = Arc::clone(&transposons_arc);
        pool.execute(move || {
            loop {
                let mut te_alignment_read = String::new();
                // read a line into te_alignment_read, then
                // break if reached EOF, else do nothing

                // block if another reader has the mutex
                let mut te_aligned_reader_child = te_aligned_reader_arc_clone.lock().unwrap();
                match te_aligned_reader_child.read_line(&mut te_alignment_read) {
                    Err(_) => panic!("Something went wrong - unable to read file"),
                    Ok(0) => break,
                    Ok(_) => (),
                }
                // unlock the mutex
                std::mem::drop(te_aligned_reader_child);
                // update the line number
                let mut line_num_child = line_num_arc_clone.lock().unwrap();
                *line_num_child += 1;
                // print status every 100,000 lines
                if *line_num_child % 100_000 == 0 {
                    println!("processing line: {}", line_num_child);
                }
                // unlock the mutex
                std::mem::drop(line_num_child);

                // create returns None if we get an unmapped read
                // or a non-split read
                // Some is only returned for split reads
                // TeAlignment's are automatically formatted in fasta format
                // with all the required info
                match TeAlignment::create(te_alignment_read, &transposons_arc_clone) {
                    None => {
                        continue;
                    }
                    Some(alignment) => {
                        let mut selected_reads_writer_child =
                            selected_reads_writer_arc_clone.lock().unwrap();
                        selected_reads_writer_child
                            .write(format!("{}\n", alignment).as_bytes())
                            .unwrap();
                        // unlock the mutex
                        std::mem::drop(selected_reads_writer_child);
                    }
                }
            }
        });
    }
    pool.join();
    return transposons_clone;
}
