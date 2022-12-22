use path_abs::PathFile;
use serde_json;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

use super::genome_alignment::GenomeAlignment;
use super::output_data_types::OutputInsertions;
use super::second_sam_file;

pub fn select_alignments(
    chroms: Vec<String>,
    min_tsd_length: u64,
    max_tsd_length: u64,
    min_te_length: f64,
    max_te_length: f64,
    genome_aligned_path: &PathFile,
    output_path: &PathFile,
    transposons_map: &HashMap<String, u64>,
    output_should_be_json: bool,
) {
    let mut second_sam_file_reader = BufReader::new(File::open(genome_aligned_path).unwrap());
    let mut output_writer = BufWriter::new(File::create(output_path).unwrap());
    second_sam_file::skip_all_comments(&mut second_sam_file_reader);
    let mut bin_heaps =
        second_sam_file::read_all_alignments_into_bin_heaps(&mut second_sam_file_reader, &chroms);
    let mut output: Vec<OutputInsertions> = Vec::new();
    for chrom in chroms {
        let non_reference = GenomeAlignment::get_non_ref_tes(
            &mut bin_heaps.get_mut(&chrom).unwrap().0,
            min_tsd_length,
            max_tsd_length,
            &chrom,
        );
        let reference = GenomeAlignment::get_ref_tes(
            &mut bin_heaps.get_mut(&chrom).unwrap().1,
            min_te_length,
            max_te_length,
            &transposons_map,
            &chrom,
        );
        output.push(OutputInsertions {
            non_reference,
            reference,
        });
    }

    if output_should_be_json {
        output_writer
            .write_all(serde_json::to_string_pretty(&output).unwrap().as_bytes())
            .unwrap();
    } else {
        output_writer.write("Chromosome\tTSD Upstream\tTSD Downstream\tOrientation\tName\t# Upstream Reads\t# Downstream Reads\tFound in Reference?\n".as_bytes()).unwrap();
        for chrom in output {
            for insertion in chrom.non_reference {
                output_writer
                    .write_all(format!("{}\n", insertion).as_bytes())
                    .unwrap();
            }
            for insertion in chrom.reference {
                output_writer
                    .write_all(format!("{}\n", insertion).as_bytes())
                    .unwrap();
            }
        }
    }
}
