// process a split-read using Regex

// CIGAR string [numbers]letter[numbers]letter...
// Each set of numbers modifies the letter after it
// S means "no match", M means "match"
// We are only interested in the 2 letters S and M
// note that H also means "no match" but it uses hard clipping
// Matches with hard clipping are only useful to us in phase 3 (process alignments)
// because they are not confident / could match in multiple spots.

// In select_reads (first selection step),
// we are interested in CIGAR strings that are
// \d+S\d+M$
// (a non-match followed by a match)
// and
// \d+M\d+S$
// (a match followed by a non-match)

// In process_reads (second selection step),
// we are also interested in these CIGAR strings

// If the string starts with a match, then
// POS is the position of the first nucleotide
// If the string does not start with a match, then
// POS is the position of the first MATCHING nucleotide

// If the match is +/+ or +/-, the CIGAR string ALWAYS starts
// from the 5' end and goes towards the 3' end

// In other words, if the match is +/- then BWA just flips it
// and continues pretending that the flipped read is the actual thing

// It turns out that by a happy accident,
// if the CIGAR string in the TE alignment step is _S_M
// and the CIGAR string in the genome alignment is _M_S,
// or vice versa, then the transposon is inserted +/+

// If the CIGAR string in the both alignments is _M_S
// or the CIGAR string in both alignments is _S_M,
// then the transposon is inserted +/-

// For this reason, it doesn't actually matter to us
// if the alignment itself is +/+ or +/-,
// it just matters if it is SM or MS,
// and which end of the transposon it aligned to (start or end)
// where "start" is the 5' end and "end" is the 3' end

use lazy_static::lazy_static;
use regex::Regex;

// lazy static regexes to avoid multiple compilations
lazy_static! {
    static ref SM_REGEX: Regex = Regex::new(r"^(\d+)S(\d+)M$").unwrap();
    static ref MS_REGEX: Regex = Regex::new(r"^(\d+)M(\d+)S$").unwrap();
    static ref HM_REGEX: Regex = Regex::new(r"^(\d+)H(\d+)M$").unwrap();
    static ref MH_REGEX: Regex = Regex::new(r"^(\d+)M(\d+)H$").unwrap();
    static ref M_REGEX: Regex = Regex::new(r"^\d+M$").unwrap();
}

#[derive(Debug)]
pub struct SMAlignment {
    pub s: u64,
    pub m: u64,
    pos: u64,
}

impl SMAlignment {
    /* never used
    pub fn get_last_s(&self) -> u64 {
        return self.pos - 1;
    }
    */
    pub fn get_first_m(&self) -> u64 {
        return self.pos;
    }
}

#[derive(Debug)]
pub struct MSAlignment {
    pub m: u64,
    pub s: u64,
    pos: u64,
}

impl MSAlignment {
    pub fn get_last_m(&self) -> u64 {
        return self.pos + self.m - 1;
    }
    /* never used
    pub fn get_first_s(&self) -> u64 {
        return self.pos + self.m;
    }
    */
}

#[derive(Debug)]
pub struct MAlignment {
    // old_s and old_m are from the TE alignment
    pub old_s: u64,
    pub old_m: u64,
    // is the TE alignment at the start of the transposon, or the end?
    is_start: bool,
    // is the genome alignment +/+? tells us the orientation of the insertion
    new_plus: bool,
    // the position of the genome alignment
    new_pos: u64,
}

impl MAlignment {
    // the old M on the boundary is either the first or last nucleotide of the transposon
    pub fn get_boundary_old_m(&self) -> u64 {
        if self.new_plus {
            // start => the TE match is SM
            if self.is_start {
                return self.new_pos + self.old_s;
            }
            // end => the TE match is MS
            else {
                return self.new_pos + self.old_m - 1;
            }
        } else {
            if self.is_start {
                return self.new_pos + self.old_m - 1;
            } else {
                return self.new_pos + self.old_s;
            }
        }
    }
}

#[derive(Debug)]
pub enum SplitReadTE {
    SM(SMAlignment),
    MS(MSAlignment),
}

// in SplitReadGenome, we parse H as if it were S
#[derive(Debug)]
pub enum SplitReadGenome {
    SM(SMAlignment),
    MS(MSAlignment),
    M(MAlignment),
}

// get a numeric capture from a regex less verbosely
// steps:
// 1. unwrap the capture (assume that it exists)
// 2. get the Match object of the capture in position "which"
// 3. unwrap (assume that it exists, since the Regex has already been checked)
// 4. turn the Match into an &str
// 5. parse out the u64
// 6. unwrap (assume that it is a number, which we know because of the Regex itself)
fn get_capture<'t>(captures: Option<regex::Captures<'t>>, which: usize) -> u64 {
    captures
        .unwrap()
        .get(which)
        .unwrap()
        .as_str()
        .parse()
        .unwrap()
}

impl SplitReadTE {
    pub fn parse(cigar: String, pos: u64) -> Option<SplitReadTE> {
        if SM_REGEX.is_match(&cigar[..]) {
            let s: u64 = get_capture(SM_REGEX.captures(&cigar[..]), 1);
            let m: u64 = get_capture(SM_REGEX.captures(&cigar[..]), 2);
            Some(SplitReadTE::SM(SMAlignment {
                s: s,
                m: m,
                pos: pos,
            }))
        } else if MS_REGEX.is_match(&cigar[..]) {
            let m: u64 = get_capture(MS_REGEX.captures(&cigar[..]), 1);
            let s: u64 = get_capture(MS_REGEX.captures(&cigar[..]), 2);
            Some(SplitReadTE::MS(MSAlignment {
                m: m,
                s: s,
                pos: pos,
            }))
        } else {
            None
        }
    }
}

impl SplitReadGenome {
    pub fn parse(
        cigar: String,
        old_s: u64,
        old_m: u64,
        is_start: bool,
        new_plus: bool,
        pos: u64,
    ) -> Option<SplitReadGenome> {
        if HM_REGEX.is_match(&cigar[..]) {
            let h: u64 = get_capture(HM_REGEX.captures(&cigar[..]), 1);
            let m: u64 = get_capture(HM_REGEX.captures(&cigar[..]), 2);
            Some(SplitReadGenome::SM(SMAlignment {
                s: h,
                m: m,
                pos: pos,
            }))
        } else if MH_REGEX.is_match(&cigar[..]) {
            let m: u64 = get_capture(MH_REGEX.captures(&cigar[..]), 1);
            let h: u64 = get_capture(MH_REGEX.captures(&cigar[..]), 2);
            Some(SplitReadGenome::MS(MSAlignment {
                m: m,
                s: h,
                pos: pos,
            }))
        } else if SM_REGEX.is_match(&cigar[..]) {
            let s: u64 = get_capture(SM_REGEX.captures(&cigar[..]), 1);
            let m: u64 = get_capture(SM_REGEX.captures(&cigar[..]), 2);
            Some(SplitReadGenome::SM(SMAlignment {
                s: s,
                m: m,
                pos: pos,
            }))
        } else if MS_REGEX.is_match(&cigar[..]) {
            let m: u64 = get_capture(MS_REGEX.captures(&cigar[..]), 1);
            let s: u64 = get_capture(MS_REGEX.captures(&cigar[..]), 2);
            Some(SplitReadGenome::MS(MSAlignment {
                m: m,
                s: s,
                pos: pos,
            }))
        } else if M_REGEX.is_match(&cigar[..]) {
            Some(SplitReadGenome::M(MAlignment {
                old_s: old_s,
                old_m: old_m,
                is_start: is_start,
                new_plus: new_plus,
                new_pos: pos,
            }))
        } else {
            None
        }
    }
}