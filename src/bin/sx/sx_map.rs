use path_abs::{PathDir, PathFile, PathOps};

use crate::te_mapper_utils::{select_alignments, select_reads};
use crate::utils;
use crate::utils::Reads;

pub fn map(
    ref_name: &str,
    reads: &Reads,
    transposons_name: &str,
    result_dir: &str,
    bwa_threads: u16,
    output_should_be_json: bool,
) {
    // create the result directory if it's not already there
    match PathDir::create(result_dir) {
        Ok(_) => (),
        Err(e) => panic!("Unable to create result directory: {}", e),
    };
    // make sure that the reference and transposons files are present
    match PathFile::new(ref_name) {
        Ok(_) => (),
        Err(e) => panic!("Reference genome not present: {}", e),
    };
    match PathFile::new(transposons_name) {
        Ok(_) => (),
        Err(e) => panic!("Transposons file not present: {}", e),
    };
    // index the transposons file and reference sequence if necessary
    utils::bwa_index_if_required(transposons_name);
    utils::bwa_index_if_required(ref_name);
    // phase 1: align the reads to the transposons
    println!("\n\nPHASE 1\n");
    let te_aligned_name = format!("{}/te_aligned.sam", result_dir);
    utils::bwa_mem_align(transposons_name, reads, &te_aligned_name[..], bwa_threads);
    // phase 2: look for split-reads (reads that go off one end of the transposon)
    // in order to be safe, only perfect matches are used
    println!("\n\nPHASE 2\n");
    let te_aligned_path = PathFile::new(te_aligned_name).unwrap();
    let result_dir_path = PathDir::new(result_dir).unwrap();
    let selected_reads_path =
        PathFile::create(result_dir_path.concat("selected_reads.fasta").unwrap()).unwrap();
    let transposons_map = select_reads::select_reads(&te_aligned_path, &selected_reads_path);
    // phase 3: align the potential split-reads to the genome and make sure that
    // the other half of the split-read is a perfect match as well
    println!("\n\nPHASE 3\n");
    let selected_reads_name = format!("{}/selected_reads.fasta", result_dir);
    let genome_aligned_name = format!("{}/genome_aligned.sam", result_dir);
    utils::bwa_mem_align(
        ref_name,
        &Reads::SingleEnd(selected_reads_name),
        &genome_aligned_name[..],
        bwa_threads,
    );
    // phase 4: select the alignments that are properly positioned on a break-point
    // between a transposon and the genome (down to the exact nucleotide)
    println!("\n\nPHASE 4\n");
    let genome_aligned_path = PathFile::new(genome_aligned_name).unwrap();

    let output_path;
    if output_should_be_json {
        output_path =
            PathFile::create(result_dir_path.concat("te_mapper_output.json").unwrap()).unwrap();
    } else {
        output_path =
            PathFile::create(result_dir_path.concat("te_mapper_output.tsv").unwrap()).unwrap();
    }

    // Drosophila Melanogaster has these 7 chromosomes (change them for a different organism)
    let chroms = vec![
        "2L".to_owned(),
        "2R".to_owned(),
        "3L".to_owned(),
        "3R".to_owned(),
        "4".to_owned(),
        "X".to_owned(),
        "Y".to_owned(),
    ];

    // params (you can change these depending on the situation)
    // min TSD length: 0
    // max TSD length: 100
    // min TE length (for reference TE's): 0.1 * the original length
    // max TE length (for reference TE's): 1.5 * the original length
    select_alignments::select_alignments(
        chroms,
        0,
        100,
        0.1,
        1.5,
        &genome_aligned_path,
        &output_path,
        &transposons_map,
        output_should_be_json,
    );
    println!("\n\nTE mapping done\n");
}
