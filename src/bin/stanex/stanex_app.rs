// the CLI app and sub commands

use clap::{App, AppSettings, Arg, SubCommand};

// the download subcommand
fn download_sc() -> App<'static, 'static> {
    SubCommand::with_name("download")
        .about("Download any file (e.g. reference sequence)")
        .arg(
            Arg::with_name("URL")
                .short("u")
                .long("url")
                .takes_value(true)
                .value_name("URL")
                .help("the URL of the file to be downloaded (e.g. ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.30_FB2019_05/fasta/dmel-all-chromosome-r6.30.fasta.gz for FlyBase release 6.30 used in Library_Prep_2020)")
                .required(true),
        )
        .arg(
            Arg::with_name("Output File")
                .short("o")
                .long("output")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to your output file (relative or absolute)")
                .required(true)
        )
}

// the variants subcommand
fn variants_sc() -> App<'static, 'static> {
    SubCommand::with_name("variants")
        .about("Find variants within a FASTQ file using BWA MEM, SAMTools, and Freebayes")
        .arg(
            Arg::with_name("Reference")
                .long("ref")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the reference sequence FASTA file (relative or absolute)")
                .required(true),
        )
        .arg(
            Arg::with_name("Paired-Ends")
                .long("paired")
                .takes_value(false)
                .help("use this argument if you want to run the variant finder in paired-ends mode")
                .required(false)
                .requires_all(&["Reads1", "Reads2"]),
                
        )
        .arg(
            Arg::with_name("Reads")
                .long("reads")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the single-end reads FASTQ file (relative or absolute)")
                .conflicts_with("Paired-Ends"),
        )
        .arg(
            Arg::with_name("Reads1")
                .long("reads1")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the 1st of 2 reads FASTQ files for paired-ends (relative or absolute)")
        )
        .arg(
            Arg::with_name("Reads2")
                .long("reads2")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the 2nd of 2 reads FASTQ files for paired-ends (relative or absolute)")
        )
        .arg(
            Arg::with_name("Result Directory")
            .long("result")
            .takes_value(true)
            .value_name("DIR")
            .help("the path to the directory where results (such as variants, alignments, and average alignment depth) will be stored (relative or absolute)")
            .required(true),
        )
        .arg(
            Arg::with_name("BWA Threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .value_name("NUM_THREADS")
                .help("the number of threads to run BWA with (default value 1; choose 1 if you want a deterministic output; choose higher numbers to run faster while taking up more memory)")
                .required(false),
        )
}

// the TE mapper subcommand
fn mapper_sc() -> App<'static, 'static> {
    SubCommand::with_name("map")
        .about("Map TE's within a genome and a reference sequence")
        .arg(
            Arg::with_name("Reference")
                .long("ref")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the reference sequence FASTA file (relative or absolute)")
                .required(true),
        )
        .arg(
            Arg::with_name("Paired-Ends")
                .long("paired")
                .takes_value(false)
                .help("use this argument if you want to run the variant finder in paired-ends mode")
                .required(false)
                .requires_all(&["Reads1", "Reads2"]),
                
        )
        .arg(
            Arg::with_name("JSON")
                .short("j")
                .long("json")
                .takes_value(false)
                .help("use this argument if you want results to be printed in JSON (useful when passing output as input to other programs, or just for convenience)")
                .required(false),
        )
        .arg(
            Arg::with_name("Reads")
                .long("reads")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the single-end reads FASTQ file (relative or absolute)")
                .conflicts_with("Paired-Ends"),
        )
        .arg(
            Arg::with_name("Reads1")
                .long("reads1")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the 1st of 2 reads FASTQ files for paired-ends (relative or absolute)")
        )
        .arg(
            Arg::with_name("Reads2")
                .long("reads2")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the 2nd of 2 reads FASTQ files for paired-ends (relative or absolute)")
        )
        .arg(
            Arg::with_name("Transposons File")
                .long("transposons")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the file containing the transposons that the TE mapper should look for (relative or absolute)")
                .required(true)
        )
        .arg(
            Arg::with_name("Result Directory")
            .long("result")
            .takes_value(true)
            .value_name("DIR")
            .help("the path to the directory where results (a TSV file containing the found transposons as well as some intermediate files) will be stored (relative or absolute)")
            .required(true),
        )
        .arg(
            Arg::with_name("BWA Threads")
                .long("bwa-threads")
                .takes_value(true)
                .value_name("NUM_THREADS")
                .help("the number of threads to run BWA with (default value 8; choose 1 if you want a deterministic output; choose higher numbers to run faster while taking up more memory)")
                .required(false),
        )
        .arg(
            Arg::with_name("TE Mapper Threads")
                .long("mapper-threads")
                .takes_value(true)
                .value_name("NUM_THREADS")
                .help("the number of threads to run the TE mapper with (default value 8; choose 1 if you want intermediate files to be deterministic; choose higher numbers to run faster while taking up more memory; output file is deterministic either way)")
                .required(false),
        )
}

// the sg (synthetic genome) subcommand
fn sg_sc() -> App<'static, 'static> {
    SubCommand::with_name("sg")
        .about("Generate a synthetic reference genome without any transposons by removing ones found by the TE mapper")
        .arg(
            Arg::with_name("Reference")
                .long("ref")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the reference sequence FASTA file (relative or absolute)")
                .required(true),
        )
        .arg(
            Arg::with_name("Transposons File")
                .long("transposons")
                .takes_value(true)
                .value_name("FILE")
                .help("the path to the file containing the transposons found by the TE mapper")
        )
        .arg(
            Arg::with_name("Result Directory")
            .long("result")
            .takes_value(true)
            .value_name("DIR")
            .help("the path to the directory where results (a TSV file containing the found transposons as well as some intermediate files) will be stored (relative or absolute)")
            .required(true),
        )
}

// the entire CLI app
pub fn app() -> App<'static, 'static> {
    App::new("StanEx Tools")
        .author("Sanath Govindarajan")
        .version("0.1")
        .about(
            "Miscellaneous tools used for Whole-Genome Sequencing analysis in the StanEx project",
        )
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommands(vec![download_sc(), variants_sc(), mapper_sc(), sg_sc()])
}
