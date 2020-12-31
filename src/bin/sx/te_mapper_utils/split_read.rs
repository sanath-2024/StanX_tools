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

#[derive(Debug)]
pub struct SMAlignment {
    pub s: u64,
    pub m: u64,
    pub pos: u64,
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
    pub pos: u64,
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
    pub is_start: bool,
    // is the genome alignment +/+? tells us the orientation of the insertion
    pub new_plus: bool,
    // the position of the genome alignment
    pub new_pos: u64,
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
