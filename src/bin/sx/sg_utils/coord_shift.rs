use super::iloc::ILoc;

// utilities to shift coordinates from normal coords to
// TE-removed coords and vice versa

// the coords may not be translatable
// from the reference genome into the synthetic genome
// because they could be within a transposon
// If within a transposon, we want to store
// the start and end positions of that transposon
// within the reference genome

enum SGCoords {
    OutsideTransposon(u64),
    WithinTransposon(u64, u64),
}

pub struct UniversalCoords {
    chrom: String,
    normal_ref: u64,
    sg_ref: SGCoords,
}

impl UniversalCoords {
    pub fn new_from_normal(
        chrom: String,
        normal_pos: u64,
        transposons: &Vec<ILoc>,
    ) -> UniversalCoords {
        // we know that the ILoc's are sorted by chromosome and
        // are in ascending order
        // so we can also look through them in order
        let mut total_to_subtract: u64 = 0;
        for iloc in transposons {
            if iloc.chrom == chrom {
                // transposon is upstream of the nucleotide in question
                if iloc.downstream_pos < normal_pos {
                    total_to_subtract += iloc.length();
                }
                // nucleotide in question is within the transposon:
                // return "within transposon" coordinates
                else if iloc.upstream_pos < normal_pos {
                    return UniversalCoords {
                        chrom: chrom,
                        normal_ref: normal_pos,
                        sg_ref: SGCoords::WithinTransposon(iloc.upstream_pos, iloc.downstream_pos),
                    };
                }
                // transposon is downstream: break from the loop
                else {
                    break;
                }
            }
        }
        // nucleotide in question is not within a transposon:
        // return good coordinates
        UniversalCoords {
            chrom: chrom,
            normal_ref: normal_pos,
            sg_ref: SGCoords::OutsideTransposon(normal_pos - total_to_subtract),
        }
    }
    pub fn new_from_sg(chrom: String, sg_pos: u64, transposons: &Vec<ILoc>) -> UniversalCoords {
        let mut normal_nt_pos: u64 = sg_pos;
        for iloc in transposons {
            if iloc.chrom == chrom {
                // if the transposon is upstream, add it to the normal position
                // note: unlike new_from_normal, in this case, we have to constantly
                // update the value that we are checking against (since we have
                // to add in the coordinates of all the transposons to the position
                // that we check against as well)
                if iloc.downstream_pos < normal_nt_pos {
                    normal_nt_pos += iloc.length();
                }
                // if the transposon is downstream, return
                // (since transposon insertions are already sorted)
                else {
                    return UniversalCoords {
                        chrom: chrom,
                        normal_ref: normal_nt_pos,
                        sg_ref: SGCoords::OutsideTransposon(sg_pos),
                    };
                }
            }
        }
        // there are no transposons downstream of the position
        UniversalCoords {
            chrom: chrom,
            normal_ref: normal_nt_pos,
            sg_ref: SGCoords::OutsideTransposon(sg_pos),
        }
    }
}
