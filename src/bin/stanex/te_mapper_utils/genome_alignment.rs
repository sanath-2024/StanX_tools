use anyhow::{bail, Result};

use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap};

use super::output_data_types::{NonRefTE, Orientation, RefTE};
use crate::tabular::Data;

// module with some helper structs and functions to represent split reads
mod split_read_genome {
    use anyhow::{bail, Result};

    use super::super::split_read::{MAlignment, MSAlignment, SMAlignment};
    use crate::regexes;

    #[derive(Debug)]
    pub enum SplitReadGenome {
        SM(SMAlignment),
        MS(MSAlignment),
        M(MAlignment),
    }

    // in SplitReadGenome, we parse H as if it were S
    impl SplitReadGenome {
        pub fn parse(
            cigar: String,
            old_m: u64,
            old_s: u64,
            is_start: bool,
            new_plus: bool,
            pos: u64,
        ) -> Result<SplitReadGenome> {
            if regexes::HM_REGEX.is_match(&cigar[..]) {
                let h: u64 = regexes::get_capture(regexes::HM_REGEX.captures(&cigar[..]), 1);
                let m: u64 = regexes::get_capture(regexes::HM_REGEX.captures(&cigar[..]), 2);
                Ok(SplitReadGenome::SM(SMAlignment {
                    s: h,
                    m: m,
                    pos: pos,
                }))
            } else if regexes::MH_REGEX.is_match(&cigar[..]) {
                let m: u64 = regexes::get_capture(regexes::MH_REGEX.captures(&cigar[..]), 1);
                let h: u64 = regexes::get_capture(regexes::MH_REGEX.captures(&cigar[..]), 2);
                Ok(SplitReadGenome::MS(MSAlignment {
                    m: m,
                    s: h,
                    pos: pos,
                }))
            } else if regexes::SM_REGEX.is_match(&cigar[..]) {
                let s: u64 = regexes::get_capture(regexes::SM_REGEX.captures(&cigar[..]), 1);
                let m: u64 = regexes::get_capture(regexes::SM_REGEX.captures(&cigar[..]), 2);
                Ok(SplitReadGenome::SM(SMAlignment {
                    s: s,
                    m: m,
                    pos: pos,
                }))
            } else if regexes::MS_REGEX.is_match(&cigar[..]) {
                let m: u64 = regexes::get_capture(regexes::MS_REGEX.captures(&cigar[..]), 1);
                let s: u64 = regexes::get_capture(regexes::MS_REGEX.captures(&cigar[..]), 2);
                Ok(SplitReadGenome::MS(MSAlignment {
                    m: m,
                    s: s,
                    pos: pos,
                }))
            } else if regexes::M_REGEX.is_match(&cigar[..]) {
                Ok(SplitReadGenome::M(MAlignment {
                    old_s: old_s,
                    old_m: old_m,
                    is_start: is_start,
                    new_plus: new_plus,
                    new_pos: pos,
                }))
            } else {
                bail!("CIGAR string is not HM, MH, SM, MS, or M");
            }
        }
    }
}

pub use split_read_genome::SplitReadGenome;

// store all relevant info from a genome alignment
// (including the previous info from the TE alignment)
#[derive(Debug)]
pub struct GenomeAlignment {
    te_name: String,
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
    fn is_mapped(flag: String) -> bool {
        let sam_flag: u16 = flag.parse().unwrap();
        (sam_flag & 4) == 0
    }
    // is the alignment +/+? The 5th least-significant bit of the SAM flag
    fn is_plus(flag: String) -> bool {
        let sam_flag: u16 = flag.parse().unwrap();
        (sam_flag & 16) == 0
    }
    // does the alignment occur on a chromosome that we care about?
    pub fn validate_chrom(chrom: &String, chroms: &Vec<String>) -> bool {
        chroms.iter().find(|x| chrom == *x) != None
    }

    // create a TE alignment from a string (skip if it doesn't meet criteria)
    // we have 3 criteria:
    // 1. SAM flag does not "&" with 4 (4 means unmapped)
    // 2. The chromosome is an actual chromosome (like 2L and 2R)
    // 3. read is a split-read from the genome side (we already know it is from the transposon side)
    pub fn create(
        genome_alignment_data: Data,
        te_alignment_data: Data,
        chroms: &Vec<String>,
    ) -> Result<(String, GenomeAlignment)> {
        if !GenomeAlignment::is_mapped(genome_alignment_data.get("FLAG")?) {
            bail!("unmapped read");
        }

        let te_name = te_alignment_data.get("TE_NAME")?;
        let old_m: u64 = te_alignment_data.get("OLD_M")?.parse()?;
        let old_s: u64 = te_alignment_data.get("OLD_S")?.parse()?;
        let is_sm_te = te_alignment_data.get("OLD_SM")? == "SM";
        let is_start = te_alignment_data.get("START_OF_TE")? == "start";

        if is_sm_te != is_start {
            let sm_str = {
                if is_sm_te {
                    "SM"
                } else {
                    "MS"
                }
            };
            let start_str = {
                if is_start {
                    "start"
                } else {
                    "end"
                }
            };
            panic!(format!(
                "TE mapper error: TE alignment was {} and aligned to the {} of the transposon",
                sm_str, start_str
            ));
        }

        let flag = genome_alignment_data.get("FLAG")?;
        let chrom = genome_alignment_data.get("RNAME")?;
        let pos: u64 = genome_alignment_data.get("POS")?.parse()?;
        let cigar_str = genome_alignment_data.get("CIGAR")?;

        if !GenomeAlignment::validate_chrom(&chrom, chroms) {
            bail!("invalid chromosome");
        }

        let is_plus = GenomeAlignment::is_plus(flag);

        let split_read = SplitReadGenome::parse(cigar_str, old_m, old_s, is_start, is_plus, pos)?;

        Ok((
            chrom.clone(),
            GenomeAlignment {
                te_name: te_name,
                old_m: old_m,
                old_s: old_s,
                is_sm_te: is_sm_te,
                is_start: is_start,
                new_plus: is_plus,
                chrom: chrom,
                split_read_genome: split_read,
            },
        ))
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
        self.te_name == other.te_name && self.get_boundary_nt() == other.get_boundary_nt()
    }
}
impl Eq for GenomeAlignment {}
impl PartialOrd for GenomeAlignment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// we want the heap (a max-heap) compare from least to greatest, not greatest to least
// so reverse it while comparing
impl Ord for GenomeAlignment {
    fn cmp(&self, other: &Self) -> Ordering {
        Reverse((&self.te_name, self.get_boundary_nt()))
            .cmp(&Reverse((&other.te_name, other.get_boundary_nt())))
    }
}

impl GenomeAlignment {
    // pull the genome alignments from a binary heap and store them in a 3D vector in sorted order
    // outer dimension: which TE is it?
    // middle dimension: which position is it? (note: all insertions in the heap are in the same chromosome)
    // inner dimension: all of the reads at that specific TE insertion
    fn make_3d_vector(heap: &mut BinaryHeap<GenomeAlignment>) -> Vec<Vec<Vec<GenomeAlignment>>> {
        let mut result: Vec<Vec<Vec<GenomeAlignment>>> = Vec::new();
        let mut cur_transposon_name = "".to_string();
        let mut cur_pos: u64 = std::u64::MAX;
        while let Some(alignment) = heap.pop() {
            let new_transposon_name = &alignment.te_name;
            let new_pos = alignment.get_boundary_nt();
            // if we are looking at the same transposon
            if new_transposon_name == &cur_transposon_name {
                // if we are looking at the same position
                if new_pos == cur_pos {
                    // add the transposon to the list of alignments at the current insertion
                    result
                        .last_mut()
                        .unwrap()
                        .last_mut()
                        .unwrap()
                        .push(alignment);
                }
                // if we are looking at a different position
                else {
                    cur_pos = new_pos;
                    result.last_mut().unwrap().push(vec![alignment]);
                }
            }
            // if we are looking at a different transposon
            else {
                cur_transposon_name = new_transposon_name.clone();
                cur_pos = new_pos;
                result.push(vec![vec![alignment]]);
            }
        }
        return result;
    }

    // get the set of non-ref TE's from a binary heap of genome alignments
    // the heap will be consumed in this function
    // this function should be run once per chromosome
    pub fn get_non_ref_tes(
        alignments: &mut BinaryHeap<GenomeAlignment>,
        min_tsd_length: u64,
        max_tsd_length: u64,
        chrom_name: &String,
    ) -> Vec<NonRefTE> {
        let alignment_vector = GenomeAlignment::make_3d_vector(alignments);
        let mut tes: Vec<NonRefTE> = Vec::new();
        // each TE will have a few split-reads downstream of it,
        // and then after that will be the upstream reads
        // this is counterintuitive but due to the TSD
        for same_transposon_name in alignment_vector {
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
                                name: alignment.te_name.clone(),
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
                                            name: alignment.te_name.clone(),
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
                                        name: alignment.te_name.clone(),
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
                                name: alignment.te_name.clone(),
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
                                            name: alignment.te_name.clone(),
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
                                        name: alignment.te_name.clone(),
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

    // get the set of ref TE's from a binary heap of genome alignments
    // the heap will be consumed in this function
    // this function should be run once per chromosome
    // this function allows for insertions and deletions within the reference transposons
    pub fn get_ref_tes(
        alignments: &mut BinaryHeap<GenomeAlignment>,
        min_te_length: f64,
        max_te_length: f64,
        all_te_lengths: &HashMap<String, u64>,
        chrom_name: &String,
    ) -> Vec<RefTE> {
        let alignment_vector = GenomeAlignment::make_3d_vector(alignments);
        let mut tes: Vec<RefTE> = Vec::new();
        // each TE will have a few split-reads upstream of it,
        // and then after that will be the downstream reads
        for same_transposon_name in alignment_vector {
            let cur_te_length = *all_te_lengths
                .get(&same_transposon_name[0][0].te_name)
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
                                name: alignment.te_name.clone(),
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
                                            name: alignment.te_name.clone(),
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
                                        name: alignment.te_name.clone(),
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
                                name: alignment.te_name.clone(),
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
                                            name: alignment.te_name.clone(),
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
                                        name: alignment.te_name.clone(),
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
