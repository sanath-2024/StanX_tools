//! # TE finding algorithm:
//! step 0: get a list of reads that are aligned to TE's and also to the same chromosome
//! step 1: split by TE name; now there should be several lists, each of which is for a single TE. For each list:
//! step 2: split this list by whether it is at the start or end of the TE, and whether the genome alignment is +/+ or +/-,
//!     so there should be 4 sub-lists now
//! step 3: sort each sub-list by the location that the boundary of the TE would be at
//! step 4: form groups from each sub-list
//!     form each group by traversing the sub-list from least to greatest
//!     and grouping together all reads whose boundary locations are at most 5 nt apart (this may be changed using the parameter `group_blur`)
//! step 5: for each group, form three consensus locations (mean, median, and mode)
//! step 6: traverse the grouped sub-lists group by group, in increasing order of location, pairing the two ends of each TE

use super::genome_alignment::GenomeAlignment;
use super::output_data_types::Orientation;

use std::collections::HashMap;
use std::rc::Rc;

#[derive(Debug)]
struct ChromList {
    chrom_name: String,
    reads: Vec<GenomeAlignment>,
}

#[derive(Debug, PartialEq, Eq)]
struct TEList {
    te_name: Rc<String>,
    reads: Vec<usize>, // indexes into the original TE list
}

#[derive(Clone, Debug)]
enum TEEnd {
    Start,
    End,
}

#[derive(Debug)]
struct SubList {
    te_name: Rc<String>,
    orientation: Orientation,
    end: TEEnd,
    reads: Vec<usize>,
}

#[derive(Debug, Clone)]
struct Group {
    te_name: Rc<String>,
    orientation: Orientation,
    end: TEEnd,
    reads: Vec<usize>,
    min: u64,
    max: u64,
    mean: f64,
    median: u64,
    mode: u64,
}

#[derive(Debug)]
struct NewNonRefTE {
    pub upstream_group: Group,
    pub downstream_group: Group,
}

#[derive(Debug)]
struct NewRefTE {
    pub upstream_group: Group,
    pub downstream_group: Group,
}

struct NewAlgoTEResults {
    plus_plus_nonref: Vec<NewNonRefTE>,
    plus_plus_ref: Vec<NewRefTE>,
    plus_minus_nonref: Vec<NewNonRefTE>,
    plus_minus_ref: Vec<NewRefTE>,
}

type NewAlgoResults = HashMap<String, NewAlgoTEResults>;

fn step1(chrom_list: &mut ChromList) -> Vec<TEList> {
    chrom_list.reads.sort_by(|a, b| a.te_name.cmp(&b.te_name));
    let mut te_lists: Vec<TEList> = Vec::new();
    let mut last_te_list = TEList {
        te_name: Rc::new(chrom_list.reads[0].te_name.clone()),
        reads: vec![0],
    };
    for (index, alignment) in chrom_list.reads.iter().enumerate() {
        if index == 0 {
            continue;
        }
        if &alignment.te_name == &*last_te_list.te_name {
            last_te_list.reads.push(index);
        } else {
            if last_te_list.te_name.len() > 0 {
                te_lists.push(last_te_list);
            }
            last_te_list = TEList {
                te_name: Rc::new(alignment.te_name.clone()),
                reads: vec![index],
            }
        }
    }
    te_lists.push(last_te_list);
    return te_lists;
}

fn step2(te_list: TEList, chrom_list: &ChromList) -> (SubList, SubList, SubList, SubList) {
    let mut plus_plus_start = SubList {
        te_name: te_list.te_name.clone(),
        orientation: Orientation::PlusPlus,
        end: TEEnd::Start,
        reads: Vec::new(),
    };
    let mut plus_plus_end = SubList {
        te_name: te_list.te_name.clone(),
        orientation: Orientation::PlusPlus,
        end: TEEnd::End,
        reads: Vec::new(),
    };
    let mut plus_minus_start = SubList {
        te_name: te_list.te_name.clone(),
        orientation: Orientation::PlusMinus,
        end: TEEnd::Start,
        reads: Vec::new(),
    };
    let mut plus_minus_end = SubList {
        te_name: te_list.te_name.clone(),
        orientation: Orientation::PlusMinus,
        end: TEEnd::End,
        reads: Vec::new(),
    };
    for te_idx in te_list.reads {
        let aln = &chrom_list.reads[te_idx];
        if aln.new_plus {
            if aln.is_start {
                plus_plus_start.reads.push(te_idx);
            } else {
                plus_plus_end.reads.push(te_idx);
            }
        } else {
            if aln.is_start {
                plus_minus_start.reads.push(te_idx);
            } else {
                plus_minus_end.reads.push(te_idx);
            }
        }
    }
    return (
        plus_plus_start,
        plus_plus_end,
        plus_minus_start,
        plus_minus_end,
    );
}

fn step3(sub_list: &mut SubList, chrom_list: &ChromList) {
    sub_list
        .reads
        .sort_by_key(|te_idx| chrom_list.reads[*te_idx].get_boundary_nt());
}

fn step4(sub_list: SubList, group_blur: u64, chrom_list: &ChromList) -> Vec<Group> {
    let mut res = vec![Group {
        te_name: sub_list.te_name.clone(),
        orientation: sub_list.orientation.clone(),
        end: sub_list.end.clone(),
        reads: vec![sub_list.reads[0]],
        min: 0,
        max: 0,
        mean: 0.0,
        median: 0,
        mode: 0,
    }];
    let last_loc = chrom_list.reads[sub_list.reads[0]].get_boundary_nt();
    let mut skipped_first_elem = false;
    for read_idx in sub_list.reads {
        if !skipped_first_elem {
            skipped_first_elem = true;
            continue;
        }
        let next_read_boundary = chrom_list.reads[read_idx].get_boundary_nt();
        if next_read_boundary - last_loc <= group_blur {
            res.last_mut().unwrap().reads.push(read_idx);
        } else {
            res.push(Group {
                te_name: sub_list.te_name.clone(),
                orientation: sub_list.orientation.clone(),
                end: sub_list.end.clone(),
                reads: vec![read_idx],
                min: 0,
                max: 0,
                mean: 0.0,
                median: 0,
                mode: 0,
            });
        }
    }
    res
}

fn step5(group: &mut Group, chrom_list: &ChromList) {
    let len = group.reads.len();
    let min = chrom_list.reads[group.reads[0]].get_boundary_nt();
    let max = chrom_list.reads[group.reads[len - 1]].get_boundary_nt();
    let median = if len % 2 == 0 {
        let left_side = chrom_list.reads[group.reads[len / 2 - 1]].get_boundary_nt();
        let right_side = chrom_list.reads[group.reads[len / 2 - 1]].get_boundary_nt();
        (left_side + right_side + 1) / 2 // note: this rounds up the .5 if necessary
    } else {
        chrom_list.reads[group.reads[len / 2]].get_boundary_nt()
    };
    let mut sum = 0.0_f64;
    let mut freq: HashMap<u64, u64> = HashMap::new();
    for read_idx in &group.reads {
        let new_read_loc = chrom_list.reads[*read_idx].get_boundary_nt();
        sum += new_read_loc as f64;
        match freq.clone().get(&new_read_loc) {
            Some(val) => {
                freq.insert(new_read_loc, *val + 1);
            }
            None => {
                freq.insert(new_read_loc, 1);
            }
        }
    }
    let mean = sum / len as f64;
    let mut max_key = 0;
    let mut max_value = 0;
    for (key, value) in &freq {
        if value > &max_value {
            max_key = *key;
            max_value = *value;
        }
    }
    group.min = min;
    group.max = max;
    group.mean = mean;
    group.median = median;
    group.mode = max_key;
}

/// max_inverted_repeat should be something small but not negligible, like 20 or 30
fn step6_plus_plus_nonref(
    start_side: &Vec<Group>,
    end_side: &Vec<Group>,
    max_inverted_repeat: u64,
) -> Vec<NewNonRefTE> {
    let mut tes = Vec::new();
    let mut end_group_idx = 0;
    for start_group in start_side {
        let start_pos = start_group.median;
        loop {
            if end_group_idx >= end_side.len() {
                break;
            }
            let end_pos = end_side[end_group_idx].median;
            if end_pos > start_pos - max_inverted_repeat {
                if end_pos < start_pos {
                    tes.push(NewNonRefTE {
                        upstream_group: start_group.clone(),
                        downstream_group: end_side[end_group_idx].clone(),
                    });
                } else {
                    break;
                }
            }
            end_group_idx += 1;
        }
    }
    return tes;
}

fn step6_plus_plus_ref(
    start_side: &Vec<Group>,
    end_side: &Vec<Group>,
    min_te_length: f64,
    max_te_length: f64,
    te_lengths: &HashMap<String, u64>,
) -> Vec<NewRefTE> {
    if start_side.len() == 0 {
        return Vec::new();
    }
    let te_name = &start_side[0].te_name;
    let te_length = *te_lengths.get(&**te_name).unwrap();
    let min_length = (min_te_length * te_length as f64) as u64;
    let max_length = (max_te_length * te_length as f64) as u64;
    let mut tes = Vec::new();
    let mut end_group_idx = 0;
    for start_group in start_side {
        let start_pos = start_group.median;
        loop {
            if end_group_idx >= end_side.len() {
                break;
            }
            let end_pos = end_side[end_group_idx].median;
            if end_pos >= start_pos + min_length {
                if end_pos <= start_pos + max_length {
                    tes.push(NewRefTE {
                        upstream_group: start_group.clone(),
                        downstream_group: end_side[end_group_idx].clone(),
                    });
                } else {
                    break;
                }
            }
            end_group_idx += 1;
        }
    }
    return tes;
}

fn step6_plus_minus_nonref(
    start_side: &Vec<Group>,
    end_side: &Vec<Group>,
    max_inverted_repeat: u64,
) -> Vec<NewNonRefTE> {
    let mut tes = Vec::new();
    let mut end_group_idx = 0;
    for start_group in start_side {
        let start_pos = start_group.median;
        loop {
            if end_group_idx >= end_side.len() {
                break;
            }
            let end_pos = end_side[end_group_idx].median;
            if end_pos > start_pos {
                if end_pos < start_pos + max_inverted_repeat {
                    tes.push(NewNonRefTE {
                        upstream_group: end_side[end_group_idx].clone(),
                        downstream_group: start_group.clone(),
                    });
                } else {
                    break;
                }
            }
            end_group_idx += 1;
        }
    }
    return tes;
}

fn step6_plus_minus_ref(
    start_side: &Vec<Group>,
    end_side: &Vec<Group>,
    min_te_length: f64,
    max_te_length: f64,
    te_lengths: &HashMap<String, u64>,
) -> Vec<NewRefTE> {
    if start_side.len() == 0 {
        return Vec::new();
    }
    let te_name = &start_side[0].te_name;
    let te_length = *te_lengths.get(&**te_name).unwrap();
    let min_length = (min_te_length * te_length as f64) as u64;
    let max_length = (max_te_length * te_length as f64) as u64;
    let mut tes = Vec::new();
    let mut end_group_idx = 0;
    for start_group in start_side {
        let start_pos = start_group.median;
        loop {
            if end_group_idx >= end_side.len() {
                break;
            }
            let end_pos = end_side[end_group_idx].median;
            if end_pos >= start_pos - max_length {
                if end_pos <= start_pos - min_length {
                    tes.push(NewRefTE {
                        upstream_group: end_side[end_group_idx].clone(),
                        downstream_group: start_group.clone(),
                    });
                } else {
                    break;
                }
            }
            end_group_idx += 1;
        }
    }
    return tes;
}

macro_rules! steps345 {
    ( $chrom_list: ident; $($sub_list: ident, $sub_list_groups: ident;)+ ) => {
        $(
            step3(&mut $sub_list, $chrom_list);
            let mut $sub_list_groups = step4($sub_list, 10, $chrom_list);
            for mut group in &mut $sub_list_groups {
                step5(&mut group, $chrom_list);
            }
        )+
    }
}

fn new_algo(chrom_list: &mut ChromList, te_lengths: &HashMap<String, u64>) -> NewAlgoResults {
    let mut res = NewAlgoResults::new();
    for te_list in step1(chrom_list) {
        let te_name = (*te_list.te_name).clone();
        let sub_lists = step2(te_list, chrom_list);
        let (mut plus_plus_start, mut plus_plus_end, mut plus_minus_start, mut plus_minus_end) =
            sub_lists;
        steps345!(
            chrom_list;
            plus_plus_start, plus_plus_start_groups;
            plus_plus_end, plus_plus_end_groups;
            plus_minus_start, plus_minus_start_groups;
            plus_minus_end, plus_minus_end_groups;
        );
        let plus_plus_nonref =
            step6_plus_plus_nonref(&plus_plus_start_groups, &plus_plus_end_groups, 30);
        let plus_plus_ref = step6_plus_plus_ref(
            &plus_plus_start_groups,
            &plus_plus_end_groups,
            0.1,
            1.5,
            te_lengths,
        );
        let plus_minus_nonref =
            step6_plus_minus_nonref(&plus_minus_start_groups, &plus_minus_end_groups, 30);
        let plus_minus_ref = step6_plus_minus_ref(
            &plus_minus_start_groups,
            &plus_minus_end_groups,
            0.1,
            1.5,
            te_lengths,
        );
        res.insert(
            te_name,
            NewAlgoTEResults {
                plus_plus_nonref: plus_plus_nonref,
                plus_plus_ref: plus_plus_ref,
                plus_minus_nonref: plus_minus_nonref,
                plus_minus_ref: plus_minus_ref,
            },
        );
    }
    return res;
}

#[cfg(test)]
mod tests {
    use super::super::genome_alignment::SplitReadGenome;
    use super::super::split_read::MAlignment;
    use super::*;

    #[test]
    fn test_step1() {
        fn make_genome_alignment(idx: usize) -> GenomeAlignment {
            let te_name = (('a' as u8 as usize + idx % 7) as u8 as char).to_string();
            return GenomeAlignment {
                te_name: te_name,
                old_m: 0,
                old_s: 0,
                is_sm_te: false,
                is_start: false,
                new_plus: false,
                chrom: "2L".to_string(),
                split_read_genome: SplitReadGenome::M(MAlignment {
                    is_start: false,
                    new_plus: false,
                    old_m: 0,
                    old_s: 0,
                    new_pos: 0,
                }),
            };
        }

        // test 1: a single element
        let mut sample_chrom_list = ChromList {
            chrom_name: "2L".to_string(),
            reads: Vec::new(),
        };
        sample_chrom_list.reads.push(make_genome_alignment(0));
        let mut te_lists = step1(&mut sample_chrom_list);
        assert_eq!(
            te_lists,
            vec![TEList {
                te_name: Rc::new("a".to_string()),
                reads: vec![0]
            }]
        );

        // test 2: 5 identical elements
        sample_chrom_list = ChromList {
            chrom_name: "2L".to_string(),
            reads: Vec::new(),
        };
        for _i in 0..5 {
            sample_chrom_list.reads.push(make_genome_alignment(0));
        }
        te_lists = step1(&mut sample_chrom_list);
        assert_eq!(
            te_lists,
            vec![TEList {
                te_name: Rc::new("a".to_string()),
                reads: vec![0, 1, 2, 3, 4]
            }]
        );

        // test 3: 5 different elements
        sample_chrom_list = ChromList {
            chrom_name: "2L".to_string(),
            reads: Vec::new(),
        };
        for i in 0..5 {
            sample_chrom_list.reads.push(make_genome_alignment(i));
        }
        te_lists = step1(&mut sample_chrom_list);
        assert_eq!(
            te_lists,
            vec![
                TEList {
                    te_name: Rc::new("a".to_string()),
                    reads: vec![0]
                },
                TEList {
                    te_name: Rc::new("b".to_string()),
                    reads: vec![1]
                },
                TEList {
                    te_name: Rc::new("c".to_string()),
                    reads: vec![2]
                },
                TEList {
                    te_name: Rc::new("d".to_string()),
                    reads: vec![3]
                },
                TEList {
                    te_name: Rc::new("e".to_string()),
                    reads: vec![4]
                }
            ]
        );

        // test 4: 7 different elements that wrap around before they are sorted
        sample_chrom_list = ChromList {
            chrom_name: "2L".to_string(),
            reads: Vec::new(),
        };
        for i in 0..100 {
            sample_chrom_list.reads.push(make_genome_alignment(i % 7));
        }
        te_lists = step1(&mut sample_chrom_list);
        assert_eq!(te_lists.len(), 7);
        sample_chrom_list
            .reads
            .sort_by(|a, b| a.te_name.cmp(&b.te_name));
        for i in 0..te_lists.len() {
            let te_name = (('a' as u8 as usize + i % 7) as u8 as char).to_string();
            let mut reads: Vec<usize> = Vec::new();
            for j in 0..sample_chrom_list.reads.len() {
                if sample_chrom_list.reads[j].te_name == te_name {
                    reads.push(j);
                }
            }
            assert_eq!(*te_lists[i].te_name, te_name);
            assert_eq!(te_lists[i].reads, reads);
        }
    }
}
