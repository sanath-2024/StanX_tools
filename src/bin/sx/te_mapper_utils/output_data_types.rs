use serde::{Deserialize, Serialize};

use std::fmt::{Display, Formatter, Result};

use super::genome_alignment::SplitReadRanges;

// I could store orientation in a bool
// but this is more readable
#[derive(Eq, PartialEq, Debug, Serialize, Deserialize, Clone)]
pub enum Orientation {
    PlusPlus,
    PlusMinus,
}

// struct to transform the TE insertion info into a TSD in a coordinate system
// (see http://bergmanlab.genetics.uga.edu/?p=36 for info about coordinate systems)
// currently, one-based fully closed and zero-based half-open are implemented
// and one-based fully closed is the default (since it is the default for BWA and BLAST)
// to use a different coordinate system, simply implement the conversion
// and use that conversion instead of the default one in
// NonRefTE::get_coords and "impl fmt for NonRefTE"
// and same for RefTE
// Note: we allow dead code here in case one or more options are not being used
#[allow(dead_code)]
#[derive(Debug)]
enum TSDCoords {
    OneBasedFullyClosed { start_pos: u64, end_pos: u64 },
    ZeroBasedHalfOpen { start_pos: u64, end_pos: u64 },
}

// struct NonRefTE keeps the TE insertion info relevant to the final TSV file
// that is not already within the genome_aligned file
// the TE is NOT found in the reference
// upstream_pos is the final nucleotide which matches the genome on the 5' end (relative to the genome) of the insertion
// downstream_pos is the first nucleotide which matches the genome on the 3' end (relative to the genome) of the insertion
// Notes: upstream_pos should be greater than downstream_pos if it's non-reference
// because of the target site duplication
// upstream_pos is the last M (match to genome) in an MS match
// downstream_pos is the first M (match to genome) in a SM match
#[derive(Debug, Serialize, Deserialize)]
pub struct NonRefTE {
    pub name: String,
    pub chrom: String,
    pub upstream_pos: u64,
    pub downstream_pos: u64,
    pub orientation: Orientation,
    pub upstream_reads: Vec<SplitReadRanges>,
    pub downstream_reads: Vec<SplitReadRanges>,
}

impl NonRefTE {
    // get which nucleotides are in the tsd from a NonRefTE struct
    fn get_coords(&self) -> TSDCoords {
        // one-based fully-closed
        return TSDCoords::OneBasedFullyClosed {
            start_pos: self.downstream_pos,
            end_pos: self.upstream_pos,
        };
        // zero-based half-open
        /*
        return TSDCoords::ZeroBasedHalfOpen {
            start_pos: self.downstream_pos - 1,
            end_pos: self.upstream_pos,
        };
        */
    }
}

// how to display a non-reference TE by default
// now we change the coordinate system if needed
impl Display for NonRefTE {
    fn fmt(&self, f: &mut Formatter) -> Result {
        let orientation_string = match &self.orientation {
            Orientation::PlusPlus => "+/+",
            Orientation::PlusMinus => "+/-",
        };
        match self.get_coords() {
            TSDCoords::OneBasedFullyClosed { start_pos, end_pos } => write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chrom,
                start_pos,
                end_pos,
                orientation_string,
                self.name,
                self.upstream_reads.len(),
                self.downstream_reads.len(),
                "non-reference",
            ),
            TSDCoords::ZeroBasedHalfOpen { start_pos, end_pos } => write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chrom,
                start_pos,
                end_pos,
                orientation_string,
                self.name,
                self.upstream_reads.len(),
                self.downstream_reads.len(),
                "non-reference",
            ),
        }
    }
}

// struct RefTE keeps the TE insertion info relevant to the final TSV file
// that is not already within the genome_aligned file
// the TE IS found in the reference
// upstream_pos is the final nucleotide which matches the genome on the 5' end (relative to the genome) of the insertion
// downstream_pos is the first nucleotide which matches the genome on the 3' end (relative to the genome) of the insertion
// Notes: upstream_pos should be less than downstream_pos if it's reference
#[derive(Debug, Serialize, Deserialize)]
pub struct RefTE {
    pub name: String,
    pub chrom: String,
    pub upstream_pos: u64,
    pub downstream_pos: u64,
    pub orientation: Orientation,
    pub upstream_reads: Vec<SplitReadRanges>,
    pub downstream_reads: Vec<SplitReadRanges>,
}

impl RefTE {
    // get which nucleotides are in the tsd from a RefTE struct
    fn get_coords(&self) -> TSDCoords {
        // one-based fully-closed
        return TSDCoords::OneBasedFullyClosed {
            start_pos: self.upstream_pos,
            end_pos: self.downstream_pos,
        };
        // zero-based half-open
        /*
        return TSDCoords::ZeroBasedHalfOpen {
            start_pos: self.upstream_pos - 1,
            end_pos: self.downstream_pos,
        };
        */
    }
}

// how to display a non-reference TE by default
// now we change the coordinate system if needed
impl Display for RefTE {
    fn fmt(&self, f: &mut Formatter) -> Result {
        let orientation_string = match &self.orientation {
            Orientation::PlusPlus => "+/+",
            Orientation::PlusMinus => "+/-",
        };
        match self.get_coords() {
            TSDCoords::OneBasedFullyClosed { start_pos, end_pos } => write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chrom,
                start_pos,
                end_pos,
                orientation_string,
                self.name,
                self.upstream_reads.len(),
                self.downstream_reads.len(),
                "reference",
            ),
            TSDCoords::ZeroBasedHalfOpen { start_pos, end_pos } => write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chrom,
                start_pos,
                end_pos,
                orientation_string,
                self.name,
                self.upstream_reads.len(),
                self.downstream_reads.len(),
                "reference",
            ),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct OutputInsertions {
    pub non_reference: Vec<NonRefTE>,
    pub reference: Vec<RefTE>,
}
