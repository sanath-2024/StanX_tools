use path_abs::PathFile;
use serde_json;

use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use super::preprocess_alignments::{GenomeAlignment, NonRefTE, RefTE};
use super::split_read::SplitReadGenome;

pub fn select_alignments(
    chroms: Vec<String>,
    min_tsd_length: u64,
    max_tsd_length: u64,
    min_te_length: f64,
    max_te_length: f64,
    genome_aligned_path: &PathFile,
    tsv_output_path: &PathFile,
    transposons_map: &HashMap<String, u64>,
    output_should_be_json: bool,
) {
    // store the split-reads in a map of heaps so that
    // they can be pulled one by one in increasing order
    let mut all_alignment_heaps = HashMap::new();
    for chrom in chroms.iter() {
        let non_ref_bin_heap: BinaryHeap<GenomeAlignment> = BinaryHeap::new();
        let ref_bin_heap: BinaryHeap<GenomeAlignment> = BinaryHeap::new();
        all_alignment_heaps.insert(chrom, vec![non_ref_bin_heap, ref_bin_heap]);
    }

    // store the split-reads from least to greatest in the heap
    for line_result in BufReader::new(File::open(genome_aligned_path).unwrap()).lines() {
        let line = line_result.unwrap();
        // check that the line is not a comment
        if line.chars().nth(0).unwrap() != '@' {
            if let Some((chrom, alignment)) = GenomeAlignment::create(line, &chroms) {
                match alignment.split_read_genome {
                    // ref
                    SplitReadGenome::M(_) => {
                        all_alignment_heaps.get_mut(&chrom).unwrap()[1].push(alignment);
                    }
                    // non-ref
                    _ => {
                        all_alignment_heaps.get_mut(&chrom).unwrap()[0].push(alignment);
                    }
                }
            }
        }
    }

    let mut final_output_writer = BufWriter::new(File::create(&tsv_output_path).unwrap());
    // first, write the TSV header if necessary
    if output_should_be_json {
        final_output_writer.write("Chromosome\tTSD Upstream\tTSD Downstream\tOrientation\tName\t# Upstream Reads\t# Downstream Reads\tFound in Reference?\n".as_bytes()).unwrap();
    } else {
        final_output_writer.write("[\n".as_bytes()).unwrap();
    }

    // then, write all the insertions
    for chrom in &chroms {
        let mut nonref_chrom_alignment_heap = &mut all_alignment_heaps.get_mut(&chrom).unwrap()[0];
        let non_ref_insertions = NonRefTE::get_non_ref_tes(
            &mut nonref_chrom_alignment_heap,
            min_tsd_length,
            max_tsd_length,
            chrom,
        );
        let mut ref_chrom_alignment_heap = &mut all_alignment_heaps.get_mut(&chrom).unwrap()[1];
        let ref_insertions = RefTE::get_ref_tes(
            &mut ref_chrom_alignment_heap,
            min_te_length,
            max_te_length,
            &transposons_map,
            chrom,
        );
        if output_should_be_json {
            for insertion in non_ref_insertions {
                final_output_writer
                    .write(
                        format!("{}\n", serde_json::to_string_pretty(&insertion).unwrap())
                            .as_bytes(),
                    )
                    .unwrap();
            }
            for insertion in ref_insertions {
                final_output_writer
                    .write(
                        format!("{}\n", serde_json::to_string_pretty(&insertion).unwrap())
                            .as_bytes(),
                    )
                    .unwrap();
            }
        } else {
            for insertion in non_ref_insertions {
                final_output_writer
                    .write(format!("{}\n", insertion).as_bytes())
                    .unwrap();
            }
            for insertion in ref_insertions {
                final_output_writer
                    .write(format!("{}\n", insertion).as_bytes())
                    .unwrap();
            }
        }
    }
}
