mod regexes;
mod sg_utils;
mod sx_app;
mod sx_download;
mod sx_map;
mod sx_variants;
mod tabular;
mod te_mapper_utils;
mod utils;

use std::error::Error;

use crate::utils::Reads;

fn main() -> Result<(), Box<dyn Error>> {
    let app = sx_app::app();
    let app_matches = app.get_matches();

    // handle "download" subcommand
    if let Some(matches) = app_matches.subcommand_matches("download") {
        let url_arg = matches.value_of("URL").unwrap();
        let output_arg = matches.value_of("Output File").unwrap();
        sx_download::download(url_arg, output_arg);
    }

    // handle "variants" subcommand
    if let Some(matches) = app_matches.subcommand_matches("variants") {
        let reference = matches.value_of("Reference").unwrap();
        let result_dir = matches.value_of("Result Directory").unwrap();
        let bwa_threads = match matches.value_of("BWA Threads") {
            Some(num) => num
                .to_owned()
                .parse::<u16>()
                .expect("Please enter a positive number of BWA threads or omit the argument"),
            None => 1,
        };
        let paired_ends = matches.is_present("Paired-Ends");
        if paired_ends {
            let reads1 = matches.value_of("Reads1").unwrap();
            let reads2 = matches.value_of("Reads2").unwrap();
            let reads_struct = Reads::PairedEnds(reads1.to_owned(), reads2.to_owned());
            sx_variants::run_variant_calling_pipeline(
                reference,
                reads_struct,
                result_dir,
                bwa_threads,
            );
        } else {
            let reads = matches.value_of("Reads").unwrap();
            let reads_struct = Reads::SingleEnd(reads.to_owned());
            sx_variants::run_variant_calling_pipeline(
                reference,
                reads_struct,
                result_dir,
                bwa_threads,
            );
        }
    }

    // handle "map" subcommand
    if let Some(matches) = app_matches.subcommand_matches("map") {
        let reference = matches.value_of("Reference").unwrap();
        let paired_ends = matches.is_present("Paired-Ends");
        let json_output = matches.is_present("JSON");
        let transposons = matches.value_of("Transposons File").unwrap();
        let result_dir = matches.value_of("Result Directory").unwrap();
        let bwa_threads = match matches.value_of("BWA Threads") {
            Some(num) => num
                .to_owned()
                .parse::<u16>()
                .expect("Please enter a positive number of BWA threads or omit the argument"),
            None => 8,
        };
        if paired_ends {
            let reads1 = matches.value_of("Reads1").unwrap();
            let reads2 = matches.value_of("Reads2").unwrap();
            let reads_struct = Reads::PairedEnds(reads1.to_owned(), reads2.to_owned());
            sx_map::map(
                reference,
                &reads_struct,
                transposons,
                result_dir,
                bwa_threads,
                json_output,
            );
        } else {
            let reads = match matches.value_of("Reads") {
                Some(reads_path) => reads_path,
                None => {
                    eprintln!("Please provide a value to the command-line argument \"reads\"");
                    std::process::exit(2);
                }
            };
            let reads_struct = Reads::SingleEnd(reads.to_owned());
            sx_map::map(
                reference,
                &reads_struct,
                transposons,
                result_dir,
                bwa_threads,
                json_output,
            );
        }
    }

    // handle "sg" subcommand
    if let Some(matches) = app_matches.subcommand_matches("sg") {
        let reference = matches.value_of("Reference").unwrap();
        // let transposons = matches.value_of("Transposons File").unwrap();
        let result_dir = matches.value_of("Result Directory").unwrap();
        sg_utils::tile_ref::tile_ref(reference, result_dir);
    }

    return Ok(());
}
