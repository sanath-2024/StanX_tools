use anyhow::Result;
use lazy_static::lazy_static;

use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

use super::genome_alignment::{GenomeAlignment, SplitReadGenome};
use crate::tabular::Metadata;

lazy_static! {
    static ref SECOND_SAM_FILE_GENOME_ALIGNMENT_METADATA: Metadata = {
        let mut headings = HashMap::new();
        headings.insert(1, "QNAME".to_string());
        headings.insert(2, "FLAG".to_string());
        headings.insert(3, "RNAME".to_string());
        headings.insert(4, "POS".to_string());
        headings.insert(6, "CIGAR".to_string());
        Metadata {
            delimiter: "\t".to_string(),
            headings: headings,
        }
    };
    static ref SECOND_SAM_FILE_TE_ALIGNMENT_METADATA: Metadata = {
        let mut headings = HashMap::new();
        headings.insert(2, "TE_NAME".to_string());
        headings.insert(3, "OLD_M".to_string());
        headings.insert(4, "OLD_S".to_string());
        headings.insert(5, "OLD_SM".to_string());
        headings.insert(6, "START_OF_TE".to_string());
        Metadata {
            delimiter: "|".to_string(),
            headings: headings,
        }
    };
}

pub fn skip_all_comments(reader: &mut BufReader<File>) {
    // skips all comments and positions the buffered reader on the first line that is an alignment
    // read the file line by line
    // get rid of comments (comments in the SAM file start with "@SQ")
    // and ignore the last comment line (starts with "@PG")
    let mut read_line;

    loop {
        read_line = String::new();
        reader.read_line(&mut read_line).unwrap();
        if read_line.chars().nth(1).unwrap() == 'P' {
            break;
        }
    }
}

pub fn read_genome_alignment(
    alignment_str: String,
    chroms: &Vec<String>,
) -> Result<(String, GenomeAlignment)> {
    let genome_alignment_data = SECOND_SAM_FILE_GENOME_ALIGNMENT_METADATA.read(alignment_str);
    let te_alignment_data =
        SECOND_SAM_FILE_TE_ALIGNMENT_METADATA.read(genome_alignment_data.get("QNAME")?);
    return GenomeAlignment::create(genome_alignment_data, te_alignment_data, chroms);
}

pub fn read_all_alignments_into_bin_heaps(
    reader: &mut BufReader<File>,
    chroms: &Vec<String>,
) -> HashMap<String, (BinaryHeap<GenomeAlignment>, BinaryHeap<GenomeAlignment>)> {
    // return a map between chromosomes and their non-ref alignments and ref alignments
    let mut sorted_result: HashMap<
        String,
        (BinaryHeap<GenomeAlignment>, BinaryHeap<GenomeAlignment>),
    > = HashMap::new();
    // read into a vector first and then convert to a binary heap
    // for O(n) performance compared to O(n log n) for inserting elements 1 by 1
    let mut unsorted_result: HashMap<String, (Vec<GenomeAlignment>, Vec<GenomeAlignment>)> =
        HashMap::new();

    for chrom in chroms {
        unsorted_result.insert(chrom.clone(), (Vec::new(), Vec::new()));
    }

    // read in all alignments into unsorted result
    let mut genome_aligned_read;

    loop {
        genome_aligned_read = String::new();
        match reader.read_line(&mut genome_aligned_read) {
            Err(_) => panic!("Something went wrong - unable to read file"),
            Ok(0) => break,
            Ok(_) => (),
        }
        if let Ok((chrom, alignment)) = read_genome_alignment(genome_aligned_read, chroms) {
            match alignment.split_read_genome {
                // ref
                SplitReadGenome::M(_) => {
                    unsorted_result.get_mut(&chrom).unwrap().1.push(alignment);
                }
                // non-ref
                _ => {
                    unsorted_result.get_mut(&chrom).unwrap().0.push(alignment);
                }
            }
        }
    }

    // sort the result
    for chrom in chroms {
        let (unsorted_nonref, unsorted_ref) = unsorted_result.remove(chrom).unwrap();
        let sorted_nonref = BinaryHeap::from(unsorted_nonref);
        let sorted_ref = BinaryHeap::from(unsorted_ref);
        sorted_result.insert(chrom.clone(), (sorted_nonref, sorted_ref));
    }

    return sorted_result;
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Write};

    use super::super::{first_sam_file, genome_alignment::GenomeAlignment};
    use super::*;

    // #[test]
    fn test_insertion_creation() {
        let mut first_sam_file_reader = BufReader::new(File::open("test/te_aligned.sam").unwrap());
        let transposons_map = first_sam_file::read_all_tes_into_map(&mut first_sam_file_reader);
        let chroms = vec![
            "2L".to_owned(),
            "2R".to_owned(),
            "3L".to_owned(),
            "3R".to_owned(),
            "4".to_owned(),
            "X".to_owned(),
            "Y".to_owned(),
        ];
        let mut second_sam_file_reader =
            BufReader::new(File::open("test/genome_aligned.sam").unwrap());
        let mut output_writer =
            BufWriter::new(File::create("test/TEST_SECOND_SAM_te_mapper_output.tsv").unwrap());
        output_writer.write("Chromosome\tTSD Upstream\tTSD Downstream\tOrientation\tName\t# Upstream Reads\t# Downstream Reads\tFound in Reference?\n".as_bytes()).unwrap();
        skip_all_comments(&mut second_sam_file_reader);
        let mut bin_heaps =
            read_all_alignments_into_bin_heaps(&mut second_sam_file_reader, &chroms);
        for chrom in chroms {
            let non_ref_insertions = GenomeAlignment::get_non_ref_tes(
                &mut bin_heaps.get_mut(&chrom).unwrap().0,
                0,
                100,
                &chrom,
            );
            let ref_insertions = GenomeAlignment::get_ref_tes(
                &mut bin_heaps.get_mut(&chrom).unwrap().1,
                0.1,
                1.5,
                &transposons_map,
                &chrom,
            );

            for insertion in non_ref_insertions {
                output_writer
                    .write(format!("{}\n", insertion).as_bytes())
                    .unwrap();
            }
            for insertion in ref_insertions {
                output_writer
                    .write(format!("{}\n", insertion).as_bytes())
                    .unwrap();
            }
        }
        // check that the test/te_mapper_output.tsv and test/TEST_SECOND_SAM_te_mapper_output.tsv are the same
        // using the "sha1sum" command in Linux
        // NOTE: this test only works if you use the old algorithm for read_all_alignments_into_bin_heaps
    }
}
