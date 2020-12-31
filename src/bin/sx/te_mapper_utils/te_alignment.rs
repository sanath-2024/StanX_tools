use anyhow::{bail, Context, Result};

use std::collections::HashMap;
use std::fmt;
use std::fmt::{Display, Formatter};

use crate::tabular::Data;

// module with some helper structs and functions to represent split reads
mod split_read_te {
    use anyhow::{bail, Result};

    use super::super::split_read::{MSAlignment, SMAlignment};
    use crate::regexes;

    #[derive(Debug)]
    pub enum SplitReadTE {
        SM(SMAlignment),
        MS(MSAlignment),
    }

    impl SplitReadTE {
        pub fn m(&self) -> u64 {
            match self {
                SplitReadTE::SM(sm_alignment) => sm_alignment.m,
                SplitReadTE::MS(ms_alignment) => ms_alignment.m,
            }
        }
        pub fn s(&self) -> u64 {
            match self {
                SplitReadTE::SM(sm_alignment) => sm_alignment.s,
                SplitReadTE::MS(ms_alignment) => ms_alignment.s,
            }
        }
        pub fn parse(cigar: String, pos: u64) -> Result<SplitReadTE> {
            if regexes::SM_REGEX.is_match(&cigar[..]) {
                let s: u64 = regexes::get_capture(regexes::SM_REGEX.captures(&cigar[..]), 1);
                let m: u64 = regexes::get_capture(regexes::SM_REGEX.captures(&cigar[..]), 2);
                Ok(SplitReadTE::SM(SMAlignment {
                    s: s,
                    m: m,
                    pos: pos,
                }))
            } else if regexes::MS_REGEX.is_match(&cigar[..]) {
                let m: u64 = regexes::get_capture(regexes::MS_REGEX.captures(&cigar[..]), 1);
                let s: u64 = regexes::get_capture(regexes::MS_REGEX.captures(&cigar[..]), 2);
                Ok(SplitReadTE::MS(MSAlignment {
                    m: m,
                    s: s,
                    pos: pos,
                }))
            } else {
                bail!("CIGAR string is not SM or MS");
            }
        }
    }
}

pub use split_read_te::SplitReadTE;

// struct TeAlignment keeps the relevant info about a TE alignment
// note: is_sm and is_start should be the same in all selected reads
#[derive(Debug)]
pub struct TeAlignment {
    pub qname: String,  // name of the read
    pub rname: String,  // name of the transposon
    pub m_size: u64,    // size of the match
    pub s_size: u64,    // size of the alignment that's outside the transposon
    pub is_sm: bool,    // is it an SM alignment (true) or an MS alignment (false)?
    pub is_start: bool, // is it at the start (true) or end (false) of the transposon?
    pub seq: String,    // the sequence of the read
}

impl TeAlignment {
    // is the alignment mapped? The 3rd least-significant bit of the SAM flag must equal 0 (1 means unmapped)
    fn is_mapped(flag: String) -> bool {
        let sam_flag: u16 = flag.parse().unwrap();
        (sam_flag & 4) == 0
    }

    // is the CIGAR string valid? In other words, is it ...S...M or ...M...S?
    fn validate_cigar_string(
        cigar_str: String,
        pos: u64,
        rname: &String,
        transposon_lengths: &HashMap<String, u64>,
    ) -> Result<SplitReadTE> {
        let split_read = SplitReadTE::parse(cigar_str, pos);

        // if it is a SM read, we need it to match at the start of the transposon
        // if it is a MS read, we need it to match at the end of the transposon
        let transposon_length = *transposon_lengths.get(rname).context(format!(
            "unable to find transposon \"{}\" in transposon list",
            rname
        ))?;

        match split_read {
            // not a split read (matches neither regex)
            Err(e) => Err(e),
            Ok(SplitReadTE::SM(sm_read)) => {
                if sm_read.get_first_m() == 1 {
                    Ok(SplitReadTE::SM(sm_read))
                } else {
                    bail!("SM read does not align to start of transposon");
                }
            }
            Ok(SplitReadTE::MS(ms_read)) => {
                if ms_read.get_last_m() == transposon_length {
                    Ok(SplitReadTE::MS(ms_read))
                } else {
                    bail!("MS read does not align to end of transposon");
                }
            }
        }
    }

    // create a TE alignment from a tabular::Data (skip if it doesn't meet criteria)
    // we have 2 criteria:
    // 1. SAM flag does not "&" with 4 (4 means unmapped)
    // 2. read aligns at the start or end of the transposon
    pub fn create(data: Data, transposon_lengths: &HashMap<String, u64>) -> Result<TeAlignment> {
        if !TeAlignment::is_mapped(data.get("FLAG")?) {
            bail!("unmapped read");
        }

        let qname = data.get("QNAME")?;
        let rname = data.get("RNAME")?;
        let pos: u64 = data.get("POS")?.parse()?;
        let cigar_str = data.get("CIGAR")?;
        let seq = data.get("SEQ")?;

        let split_read =
            TeAlignment::validate_cigar_string(cigar_str, pos, &rname, transposon_lengths)?;

        let s_size = split_read.s();
        let m_size = split_read.m();
        let is_sm: bool;
        let is_start: bool;

        match split_read {
            SplitReadTE::SM(_) => {
                is_sm = true;
                is_start = true;
            }
            SplitReadTE::MS(_) => {
                is_sm = false;
                is_start = false;
            }
        }

        Ok(TeAlignment {
            qname: qname,
            rname: rname,
            m_size: m_size,
            s_size: s_size,
            is_sm: is_sm,
            is_start: is_start,
            seq: seq,
        })
    }
}

// how to display a TE alignment by default
impl Display for TeAlignment {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
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
