use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap};
use std::fmt;

use super::split_read::SplitReadGenome;

// READING

// store all relevant info from a genome alignment
// (including the previous info from the TE alignment)
#[derive(Debug)]
pub struct GenomeAlignment {
    rname: String,
    old_m: u64,
    old_s: u64,
    is_sm_te: bool,
    is_start: bool,
    new_plus: bool,
    chrom: String,
    pub split_read_genome: SplitReadGenome,
}

impl GenomeAlignment {
    // is the alignment mapped? The 3rd least-significant bit of the SAM flag must equal 0 (1 means unmapped)
    fn is_mapped(fields: &Vec<String>) -> bool {
        let sam_flag: u16 = fields[1][..].parse().unwrap();
        (sam_flag & 4) == 0
    }
    // is the alignment +/+? The 5th least-significant bit of the SAF flag
    fn is_plus(fields: &Vec<String>) -> bool {
        let sam_flag: u16 = fields[1][..].parse().unwrap();
        (sam_flag & 16) == 0
    }
    fn parse_rname(te_aln_fields: &Vec<String>) -> String {
        te_aln_fields[1].clone()
    }
    // (old m, old s, is SM te, is start)
    fn parse_te_fields(te_aln_fields: &Vec<String>) -> (u64, u64, bool, bool) {
        let old_m: u64 = te_aln_fields[2].parse::<u64>().unwrap();
        let old_s: u64 = te_aln_fields[3].parse::<u64>().unwrap();
        (
            old_m,
            old_s,
            te_aln_fields[4] == "SM",
            te_aln_fields[5] == "start",
        )
    }
    fn parse_chrom(fields: &Vec<String>) -> String {
        fields[2].clone()
    }

    // create a TE alignment from a string (skip if it doesn't meet criteria)
    // we have 2 criteria:
    // 1. SAM flag does not "&" with 4 (4 means unmapped)
    // 2. The chromosome is an actual chromosome (like 2L and 2R)
    // 2. read is a split-read from the genome side (we already know it is from the transposon side)
    pub fn create(read_string: String, chroms: &Vec<String>) -> Option<(String, GenomeAlignment)> {
        // BWA gives us a tab-spaced file
        let fields: Vec<String> = read_string.split("\t").map(|x| x.to_owned()).collect();
        // the select-reads step produces names split by vertical bars
        let te_aln_fields: Vec<String> = fields[0].split("|").map(|x| x.to_owned()).collect();

        // criterion 1
        if !(GenomeAlignment::is_mapped(&fields)) {
            return None;
        }

        // criterion 2
        let chrom = GenomeAlignment::parse_chrom(&fields);
        if let None = chroms.iter().find(|x| x == &&chrom) {
            return None;
        }

        // criterion 3
        // use the position and CIGAR string to determine if it is a split-read
        let cigar_string = fields[5].clone();
        let pos: u64 = fields[3][..].parse().unwrap(); // position of the match
        let is_plus: bool = GenomeAlignment::is_plus(&fields); // orientation of the match
                                                               // TE alignment info
        let (old_m, old_s, is_sm_te, is_start) = GenomeAlignment::parse_te_fields(&te_aln_fields);

        let split_read_result =
            SplitReadGenome::parse(cigar_string, old_s, old_m, is_start, is_plus, pos);
        match split_read_result {
            // not a split read (matches neither regex)
            None => {
                return None;
            }
            Some(split_read) => {
                let x = Some((
                    chrom.clone(),
                    GenomeAlignment {
                        rname: GenomeAlignment::parse_rname(&te_aln_fields),
                        old_m: old_m,
                        old_s: old_s,
                        is_sm_te: is_sm_te,
                        is_start: is_start,
                        new_plus: is_plus,
                        chrom: chrom,
                        split_read_genome: split_read,
                    },
                ));
                return x;
            }
        }
    }

    // get the position of the boundary nucleotide in the genome
    // right next to the transposon
    // this must be in a one-indexed coordinate system
    // (coordinate conversions happen in a later step)
    fn get_boundary_nt(&self) -> u64 {
        match &self.split_read_genome {
            // the upstream end of the insertion (non-reference)
            SplitReadGenome::MS(alignment) => alignment.get_last_m(),
            // the downstream end of the insertion (non-reference)
            SplitReadGenome::SM(alignment) => alignment.get_first_m(),
            // reference insertion
            SplitReadGenome::M(alignment) => alignment.get_boundary_old_m(),
        }
    }

    // is this split-read upstream or downstream of the insertion site?
    // non-reference: upstream should be a MS read while downstream would be an SM read
    // reference: upstream would be if it's start and +/+, or end and +/-
    fn upstream(&self) -> bool {
        match &self.split_read_genome {
            SplitReadGenome::MS(_) => true,
            SplitReadGenome::SM(_) => false,
            SplitReadGenome::M(_) => !(self.is_start ^ self.new_plus),
        }
    }
}

// order the genome alignments first by transposon name, then by position
// for optimal grouping in a binary heap
impl PartialEq for GenomeAlignment {
    fn eq(&self, other: &Self) -> bool {
        self.rname == other.rname && self.get_boundary_nt() == other.get_boundary_nt()
    }
}
impl Eq for GenomeAlignment {}
impl PartialOrd for GenomeAlignment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// compare from least to greatest, not greatest to least
// so reverse it while comparing
impl Ord for GenomeAlignment {
    fn cmp(&self, other: &Self) -> Ordering {
        Reverse((&self.rname, self.get_boundary_nt()))
            .cmp(&Reverse((&other.rname, other.get_boundary_nt())))
    }
}

// PROCESSING / WRITING

// I could store orientation in a bool
// but this is more readable
#[derive(Eq, PartialEq, Debug)]
enum Orientation {
    PlusPlus,
    PlusMinus,
}

// struct NonRefTE keeps the TE insertion info relevant to the final BED file
// that is not already within the genome_aligned file
// the TE is NOT found in the reference
// upstream_pos is the final nucleotide which matches the genome on the 5' end (relative to the genome) of the insertion
// downstream_pos is the first nucleotide which matches the genome on the 3' end (relative to the genome) of the insertion
// Notes: upstream_pos should be greater than downstream_pos if it's non-reference
// because of the target site duplication
// upstream_pos is the last M (match to genome) in an MS match
// downstream_pos is the first M (match to genome) in a SM match
#[derive(Debug)]
pub struct NonRefTE {
    name: String,
    chrom: String,
    upstream_pos: u64,
    downstream_pos: u64,
    orientation: Orientation,
    num_upstream_reads: u64,
    num_downstream_reads: u64,
}

// struct RefTE keeps the TE insertion info relevant to the final BED file
// that is not already within the genome_aligned file
// the TE IS found in the reference
// upstream_pos is the final nucleotide which matches the genome on the 5' end (relative to the genome) of the insertion
// downstream_pos is the first nucleotide which matches the genome on the 3' end (relative to the genome) of the insertion
// Notes: upstream_pos should be less than downstream_pos if it's reference
#[derive(Debug)]
pub struct RefTE {
    name: String,
    chrom: String,
    upstream_pos: u64,
    downstream_pos: u64,
    orientation: Orientation,
    num_upstream_reads: u64,
    num_downstream_reads: u64,
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

    // get the set of non-ref TE's from a binary heap of genome alignments
    // the heap will be consumed in this function
    // this function should be run once per chromosome
    // note that everything is still one-indexed for now
    pub fn get_non_ref_tes(
        alignments: &mut BinaryHeap<GenomeAlignment>,
        min_tsd_length: u64,
        max_tsd_length: u64,
        chrom_name: &String,
    ) -> Vec<NonRefTE> {
        // group the alignments in terms of their transposon and position
        // from least to greatest for each transposon
        let mut grouped_by_te: Vec<Vec<Vec<GenomeAlignment>>> = Vec::new();
        let mut cur_transposon_name = "".to_owned();
        let mut cur_pos: u64 = std::u64::MAX;
        while let Some(alignment) = alignments.pop() {
            let new_transposon_name = alignment.rname.clone();
            let new_pos = alignment.get_boundary_nt();
            if new_transposon_name == cur_transposon_name {
                // same transposon, same position
                if new_pos == cur_pos {
                    grouped_by_te
                        .last_mut()
                        .unwrap()
                        .last_mut()
                        .unwrap()
                        .push(alignment);
                }
                // same transposon, different position
                else {
                    grouped_by_te.last_mut().unwrap().push(vec![alignment]);
                    cur_pos = new_pos;
                }
            }
            // different transposon
            else {
                grouped_by_te.push(vec![vec![alignment]]);
                cur_transposon_name = new_transposon_name;
                cur_pos = new_pos;
            }
        }
        let mut tes: Vec<NonRefTE> = Vec::new();
        // each TE will have a few split-reads downstream of it,
        // and then after that will be the upstream reads
        // this is counterintuitive but due to the TSD
        for same_transposon_name in grouped_by_te {
            for same_position in same_transposon_name {
                for alignment in same_position {
                    let position = alignment.get_boundary_nt();
                    // upstream of the transposon, MS read
                    if alignment.upstream() {
                        let orientation = if alignment.is_sm_te {
                            Orientation::PlusPlus
                        } else {
                            Orientation::PlusMinus
                        };
                        match tes.last_mut() {
                            // no TE's in the vector yet
                            None => tes.push(NonRefTE {
                                name: alignment.rname.clone(),
                                chrom: chrom_name.clone(),
                                upstream_pos: alignment.get_boundary_nt(),
                                downstream_pos: std::u64::MAX / 2,
                                orientation: orientation,
                                num_upstream_reads: 1,
                                num_downstream_reads: 0,
                            }),
                            // if there are TE's in the vector, match against the previous ones
                            Some(insertion) => {
                                if orientation == insertion.orientation {
                                    // we are still in the same insertion
                                    // if the upstream position matches
                                    // or it is between min_tsd_length and max_tsd_length after the downstream position
                                    if position == insertion.upstream_pos {
                                        insertion.num_upstream_reads += 1;
                                    } else if position >= insertion.downstream_pos + min_tsd_length
                                        && position <= insertion.downstream_pos + max_tsd_length
                                    {
                                        insertion.upstream_pos = position;
                                        insertion.num_upstream_reads += 1;
                                    }
                                    // we are in a new insertion
                                    else {
                                        tes.push(NonRefTE {
                                            name: alignment.rname.clone(),
                                            chrom: chrom_name.clone(),
                                            upstream_pos: position,
                                            downstream_pos: std::u64::MAX / 2,
                                            orientation: orientation,
                                            num_upstream_reads: 1,
                                            num_downstream_reads: 0,
                                        });
                                    }
                                }
                                // we are in a new insertion
                                else {
                                    tes.push(NonRefTE {
                                        name: alignment.rname.clone(),
                                        chrom: chrom_name.clone(),
                                        upstream_pos: position,
                                        downstream_pos: std::u64::MAX / 2,
                                        orientation: orientation,
                                        num_upstream_reads: 1,
                                        num_downstream_reads: 0,
                                    });
                                }
                            }
                        }
                    }
                    // downstream of the transposon, SM read
                    else {
                        let orientation = if alignment.is_sm_te {
                            Orientation::PlusMinus
                        } else {
                            Orientation::PlusPlus
                        };
                        match tes.last_mut() {
                            // no TE's in the vector yet
                            None => tes.push(NonRefTE {
                                name: alignment.rname.clone(),
                                chrom: chrom_name.clone(),
                                upstream_pos: std::u64::MAX / 2,
                                downstream_pos: alignment.get_boundary_nt(),
                                orientation: orientation,
                                num_upstream_reads: 0,
                                num_downstream_reads: 1,
                            }),
                            // if there are TE's in the vector, match against the previous ones
                            Some(insertion) => {
                                if orientation == insertion.orientation {
                                    // we are still in the same insertion
                                    // only if the downstream position matches
                                    if position == insertion.downstream_pos {
                                        insertion.num_downstream_reads += 1;
                                    }
                                    // we are in a new insertion
                                    else {
                                        tes.push(NonRefTE {
                                            name: alignment.rname.clone(),
                                            chrom: chrom_name.clone(),
                                            upstream_pos: std::u64::MAX / 2,
                                            downstream_pos: position,
                                            orientation: orientation,
                                            num_upstream_reads: 0,
                                            num_downstream_reads: 1,
                                        });
                                    }
                                }
                                // we are in a new insertion
                                else {
                                    tes.push(NonRefTE {
                                        name: alignment.rname.clone(),
                                        chrom: chrom_name.clone(),
                                        upstream_pos: std::u64::MAX / 2,
                                        downstream_pos: position,
                                        orientation: orientation,
                                        num_upstream_reads: 0,
                                        num_downstream_reads: 1,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
        // if the insertion does not have reads on both ends, discard it
        // (can't use iterators because of borrowing)
        let mut filtered_tes: Vec<NonRefTE> = Vec::new();
        for insertion in tes {
            if insertion.num_upstream_reads > 0 && insertion.num_downstream_reads > 0 {
                filtered_tes.push(insertion);
            }
        }
        // finally, sort by location instead of TE name
        filtered_tes.sort_by(|first, second| {
            first
                .upstream_pos
                .partial_cmp(&second.upstream_pos)
                .unwrap()
        });
        return filtered_tes;
    }
}

// how to display a non-reference TE by default
// now we change the coordinate system if needed
impl fmt::Display for NonRefTE {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
                self.num_upstream_reads,
                self.num_downstream_reads,
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
                self.num_upstream_reads,
                self.num_downstream_reads,
                "non-reference",
            ),
        }
    }
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

    // get the set of ref TE's from a binary heap of genome alignments
    // the heap will be consumed in this function
    // this function should be run once per chromosome
    // note that everything is still one-indexed for now
    // this function allows for insertions and deletions within the reference transposons
    pub fn get_ref_tes(
        alignments: &mut BinaryHeap<GenomeAlignment>,
        min_te_length: f64,
        max_te_length: f64,
        all_te_lengths: &HashMap<String, u64>,
        chrom_name: &String,
    ) -> Vec<RefTE> {
        // group the alignments in terms of their transposon and position
        // from least to greatest for each transposon
        let mut grouped_by_te: Vec<Vec<Vec<GenomeAlignment>>> = Vec::new();
        let mut cur_transposon_name = "".to_owned();
        let mut cur_pos: u64 = std::u64::MAX;
        while let Some(alignment) = alignments.pop() {
            let new_transposon_name = alignment.rname.clone();
            let new_pos = alignment.get_boundary_nt();
            if new_transposon_name == cur_transposon_name {
                // same transposon, same position
                if new_pos == cur_pos {
                    grouped_by_te
                        .last_mut()
                        .unwrap()
                        .last_mut()
                        .unwrap()
                        .push(alignment);
                }
                // same transposon, different position
                else {
                    grouped_by_te.last_mut().unwrap().push(vec![alignment]);
                    cur_pos = new_pos;
                }
            }
            // different transposon
            else {
                grouped_by_te.push(vec![vec![alignment]]);
                cur_transposon_name = new_transposon_name;
                cur_pos = new_pos;
            }
        }
        let mut tes: Vec<RefTE> = Vec::new();
        // each TE will have a few split-reads upstream of it,
        // and then after that will be the downstream reads
        for same_transposon_name in grouped_by_te {
            let cur_te_length = *all_te_lengths
                .get(&same_transposon_name[0][0].rname)
                .unwrap() as f64;
            for same_position in same_transposon_name {
                for alignment in same_position {
                    let position = alignment.get_boundary_nt();
                    let orientation = if alignment.new_plus {
                        Orientation::PlusPlus
                    } else {
                        Orientation::PlusMinus
                    };
                    // upstream of the transposon
                    if alignment.upstream() {
                        match tes.last_mut() {
                            // no TE's in the vector yet
                            None => tes.push(RefTE {
                                name: alignment.rname.clone(),
                                chrom: chrom_name.clone(),
                                upstream_pos: alignment.get_boundary_nt(),
                                downstream_pos: std::u64::MAX / 2,
                                orientation: orientation,
                                num_upstream_reads: 1,
                                num_downstream_reads: 0,
                            }),
                            // if there are TE's in the vector, match against the previous ones
                            Some(insertion) => {
                                if orientation == insertion.orientation {
                                    // we are still in the same insertion
                                    // only if the upstream position matches
                                    if position == insertion.upstream_pos {
                                        insertion.num_upstream_reads += 1;
                                    }
                                    // we are in a new insertion
                                    else {
                                        tes.push(RefTE {
                                            name: alignment.rname.clone(),
                                            chrom: chrom_name.clone(),
                                            upstream_pos: position,
                                            downstream_pos: std::u64::MAX / 2,
                                            orientation: orientation,
                                            num_upstream_reads: 1,
                                            num_downstream_reads: 0,
                                        });
                                    }
                                }
                                // we are in a new insertion
                                else {
                                    tes.push(RefTE {
                                        name: alignment.rname.clone(),
                                        chrom: chrom_name.clone(),
                                        upstream_pos: position,
                                        downstream_pos: std::u64::MAX / 2,
                                        orientation: orientation,
                                        num_upstream_reads: 1,
                                        num_downstream_reads: 0,
                                    });
                                }
                            }
                        }
                    }
                    // downstream of the transposon
                    else {
                        match tes.last_mut() {
                            // no TE's in the vector yet
                            None => tes.push(RefTE {
                                name: alignment.rname.clone(),
                                chrom: chrom_name.clone(),
                                upstream_pos: std::u64::MAX / 2,
                                downstream_pos: alignment.get_boundary_nt(),
                                orientation: orientation,
                                num_upstream_reads: 0,
                                num_downstream_reads: 1,
                            }),
                            // if there are TE's in the vector, match against the previous ones
                            Some(insertion) => {
                                if orientation == insertion.orientation {
                                    // we are still in the same insertion
                                    // if the downstream position matches
                                    // or it is between min_te_length and max_te_length after the upstream position
                                    if position == insertion.downstream_pos {
                                        insertion.num_downstream_reads += 1;
                                    } else if position
                                        >= insertion.upstream_pos
                                            + ((min_te_length * cur_te_length) as u64)
                                        && position
                                            <= insertion.upstream_pos
                                                + ((max_te_length * cur_te_length) as u64)
                                    {
                                        insertion.downstream_pos = position;
                                        insertion.num_downstream_reads += 1;
                                    }
                                    // we are in a new insertion
                                    else {
                                        tes.push(RefTE {
                                            name: alignment.rname.clone(),
                                            chrom: chrom_name.clone(),
                                            upstream_pos: std::u64::MAX / 2,
                                            downstream_pos: position,
                                            orientation: orientation,
                                            num_upstream_reads: 0,
                                            num_downstream_reads: 1,
                                        });
                                    }
                                }
                                // we are in a new insertion
                                else {
                                    tes.push(RefTE {
                                        name: alignment.rname.clone(),
                                        chrom: chrom_name.clone(),
                                        upstream_pos: std::u64::MAX / 2,
                                        downstream_pos: position,
                                        orientation: orientation,
                                        num_upstream_reads: 0,
                                        num_downstream_reads: 1,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
        // if the insertion does not have reads on both ends, discard it
        // (can't use iterators because of borrowing)
        let mut filtered_tes: Vec<RefTE> = Vec::new();
        for insertion in tes {
            if insertion.num_upstream_reads > 0 && insertion.num_downstream_reads > 0 {
                filtered_tes.push(insertion);
            }
        }
        // finally, sort by location instead of TE name
        filtered_tes.sort_by(|first, second| {
            first
                .upstream_pos
                .partial_cmp(&second.upstream_pos)
                .unwrap()
        });
        return filtered_tes;
    }
}

// how to display a non-reference TE by default
// now we change the coordinate system if needed
impl fmt::Display for RefTE {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
                self.num_upstream_reads,
                self.num_downstream_reads,
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
                self.num_upstream_reads,
                self.num_downstream_reads,
                "reference",
            ),
        }
    }
}
