use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

use super::split_read::SplitReadTE;

// find transposon name and length from first few rows of BWA output file
// note: read_string ends with a newline
pub fn create_transposon(read_string: &String, tr_map: &mut HashMap<String, u64>) {
    let fields: Vec<String> = read_string.split('\t').map(|x| x.to_owned()).collect();
    tr_map.insert(
        fields[1][3..].to_owned(),
        fields[2][3..fields[2].len() - 1].parse().unwrap(),
    );
}

// struct TeAlignment keeps the relevant info about a TE alignment
// note: is_sm and is_start should be the same in all selected reads
#[derive(Debug)]
pub struct TeAlignment {
    qname: String,  // name of the read
    rname: String,  // name of the transposon
    m_size: u64,    // size of the match
    s_size: u64,    // size of the alignment that's outside the transposon
    is_sm: bool,    // is it an SM alignment (true) or an MS alignment (false)?
    is_start: bool, // is it at the start (true) or end (false) of the transposon?
    seq: String,    // the sequence of the read
}

impl TeAlignment {
    // is the alignment mapped? The 3rd least-significant bit of the SAM flag must equal 0 (1 means unmapped)
    fn is_mapped(fields: &Vec<String>) -> bool {
        let sam_flag: u16 = fields[1][..].parse().unwrap();
        (sam_flag & 4) == 0
    }
    fn parse_qname(fields: &Vec<String>) -> String {
        fields[0].clone()
    }
    fn parse_rname(fields: &Vec<String>) -> String {
        fields[2].clone()
    }
    fn parse_seq(fields: &Vec<String>) -> String {
        fields[9].clone()
    }

    // create a TE alignment from a string (skip if it doesn't meet criteria)
    // we have 2 criteria:
    // 1. SAM flag does not "&" with 4 (4 means unmapped)
    // 2. read aligns at the start or end of the transposon
    pub fn create(
        read_string: String,
        transposons: &Arc<HashMap<String, u64>>,
    ) -> Option<TeAlignment> {
        // BWA gives us a tab-spaced file
        let fields: Vec<String> = read_string.split("\t").map(|x| x.to_owned()).collect();
        // criterion 1
        if !(TeAlignment::is_mapped(&fields)) {
            return None;
        }
        // criterion 2
        // use the position and CIGAR string to determine if it is a split-read
        let cigar_string = fields[5].clone();
        let pos: u64 = fields[3][..].parse().unwrap(); // position of the match
        let rname = TeAlignment::parse_rname(&fields); // name of the transposon
        if let Some(transposon_length) = transposons.get(&rname) {
            let split_read = SplitReadTE::parse(cigar_string, pos);
            // if it is a SM read, we need it to match at the start of the transposon
            // if it is a MS read, we need it to match at the end of the transposon
            let (s_size, m_size, is_sm, is_start) = match split_read {
                // not a split read (matches neither regex)
                None => {
                    return None;
                }
                Some(SplitReadTE::SM(sm_read)) => {
                    if sm_read.get_first_m() == 1 {
                        (sm_read.s, sm_read.m, true, true)
                    } else {
                        return None;
                    }
                }
                Some(SplitReadTE::MS(ms_read)) => {
                    if &ms_read.get_last_m() == transposon_length {
                        (ms_read.s, ms_read.m, false, false)
                    } else {
                        return None;
                    }
                }
            };

            return Some(TeAlignment {
                qname: TeAlignment::parse_qname(&fields),
                rname: rname,
                m_size: m_size,
                s_size: s_size,
                is_sm: is_sm,
                is_start: is_start,
                seq: TeAlignment::parse_seq(&fields),
            });
        } else {
            println!("{}", rname);
            panic!("unable to find transposon in transposon list");
        };
    }
}

// how to display a TE alignment by default
impl fmt::Display for TeAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let start_str = {
            if self.is_start {
                "start"
            } else {
                "end"
            }
        };
        let sm_str = {
            if self.is_sm {
                "SM"
            } else {
                "MS"
            }
        };
        write!(
            f,
            ">{}|{}|{}|{}|{}|{}\n{}",
            self.qname, self.rname, self.m_size, self.s_size, sm_str, start_str, self.seq
        )
    }
}