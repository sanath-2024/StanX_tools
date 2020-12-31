use std::process::Command;

use crate::utils;
use crate::utils::Reads;

// "fix" alignments by cleaning up read pairing information and flags
// also compress from SAM format to BAM format to save space and
// put the output in result_dir/fixed_alignments.sam
fn samtools_fixmate(result_dir: &str) {
    let input_file = format!("{}/raw_alignments.sam", result_dir);
    let output_file = format!("{}/fixed_alignments.bam", result_dir);

    println!("Waiting for samtools fixmate...");
    let mut child_proc = Command::new("samtools")
        .args(&["fixmate", "-O", "bam", &input_file[..], &output_file[..]])
        .spawn()
        .unwrap();
    let _result = child_proc.wait().unwrap();
    println!("Alignment fixing complete");
}

// sort the alignments in numerical order (Freebayes does not work unless alignments are in numerical order)
// keep everything compressed in the BAM format to save space
fn samtools_sort(result_dir: &str) {
    let input_file = format!("{}/fixed_alignments.bam", result_dir);
    let output_file = format!("{}/sorted_alignments.bam", result_dir);

    println!("Waiting for samtools sort...");
    let mut child_proc = Command::new("samtools")
        .args(&["sort", "-O", "bam", &input_file[..], "-o", &output_file[..]])
        .spawn()
        .unwrap();
    let _result = child_proc.wait().unwrap();
    println!("Alignment sorting complete");
}

// do variant calling with Freebayes
// use the --pooled-continuous flag since we are using more than 1 fly in our sample
fn freebayes_variant_call(ref_name: &str, result_dir: &str) {
    let input_file = format!("{}/sorted_alignments.bam", result_dir);
    let output_file = format!("{}/variants.vcf", result_dir);

    println!("Waiting for Freebayes...");
    let mut child_proc = Command::new("freebayes")
        .args(&[
            "--pooled-continuous",
            "--fasta-reference",
            ref_name,
            "--bam",
            &input_file[..],
            "--vcf",
            &output_file[..],
        ])
        .spawn()
        .unwrap();
    let _result = child_proc.wait().unwrap();
    println!("Variant calling complete");
}

// run the entire pipeline, one step after another
// everything must be blocking since each step depends on the previous step's output
pub fn run_variant_calling_pipeline(
    ref_name: &str,
    reads_names: Reads,
    result_dir: &str,
    bwa_threads: u16,
) {
    utils::bwa_index_if_required(ref_name);
    utils::bwa_mem_align(
        ref_name,
        &reads_names.clone(),
        &format!("{}/raw_alignments.sam", result_dir)[..],
        bwa_threads,
    );
    samtools_fixmate(result_dir);
    samtools_sort(result_dir);
    freebayes_variant_call(ref_name, result_dir);
}
