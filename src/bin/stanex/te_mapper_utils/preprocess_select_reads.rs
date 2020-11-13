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

#[cfg(test)]
mod tests {
    use super::*;

    use std::fs::File;
    use std::io::BufReader;
    use std::io::BufRead;

    use std::sync::Arc;

    #[test]
    fn test_cigar_string_parsing() {
        // test the creation of, m_size, s_size, is_sm, and is_start

        // read in the transposons
        let mut te_aligned_reader = BufReader::with_capacity(
            65_536,
            File::open("test/te_aligned.sam").unwrap(),
        );

        // read the file line by line
        let mut te_aligned_read;

        // first, get rid of comments (comments in the SAM file start with "@SQ")
        // and ignore the last comment line (starts with "@PG")
        let mut transposons: HashMap<String, u64> = HashMap::new();
        loop {
            te_aligned_read = String::new();
            te_aligned_reader
                .read_line(&mut te_aligned_read)
                .unwrap();
            if te_aligned_read.chars().nth(1).unwrap() == 'P' {
                break;
            } else {
                create_transposon(&te_aligned_read, &mut transposons);
            }
        }

        let transposons_arc = Arc::new(transposons);

        // use actual inputs from the file
        // 54S34M62S is invalid
        let input0 = "I_MADE_UP_THIS_READ	0	roo#LTR/Bel-Pao	1	60	54S34M62S	*	0	0	CCTGGCTTGGGGCGGCCGCGGGTTCGTGGCGTCGGCGCTATTTGTTCCTTGGCAGTCGGCTCTTCCTATCATTGTGAAGCAAAATTCATATGGCATTGTCTCCTAAAACTTTTCTATAGTGCCGTATTTCTATGGCGCCCACTGTGAAGN	--F-7-F7----A--F7---------7--77----J7<---77--7---A--7-7-----<A7--7F<7FAAJA7---F-<-F7<<-<----<<--<---<--F7-F-F-<JFJAF7<JFJJAJFJFFJAJJJJJJJJJJJJFJJF<AA#	NM:i:0	MD:Z:34	AS:i:34	XS:i:0";
        assert!(TeAlignment::create(String::from(input0), &transposons_arc).is_none());

        // case 1: +/+ match at start
        // flag 0 (+/+), 119S31M => the read aligned to the top strand of the transposon, and the match was at the beginning of the transposon:
        //  match ... ---->
        //              -------->
        //              <--------
        // expected: m_size = 31, s_size = 119, is_sm = true, is_start = true
        let input1 = "2L_Read_976816	0	roo#LTR/Bel-Pao	1	0	119S31M	*	0	0	ACATATGATATAAATAGCATTAAATGTTGAGTATAACGTGTCAAAGAATCCTTGGGATGAATAATAACGGAGGAAGCTGTAAATATAACCAGATTAGAAACCTATTCCTATAAACTCTCTGTTCACACATGAACACGAATATATTTAAAG	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	NM:i:0	MD:Z:31	AS:i:31	XS:i:31	XA:Z:roo#LTR/Bel-Pao,+8665,119S31M,0;";
        let result1 = TeAlignment::create(String::from(input1), &transposons_arc);
        assert!(result1.is_some());
        if let Some(res) = result1 {
            assert_eq!(res.m_size, 31);
            assert_eq!(res.s_size, 119);
            assert!(res.is_sm);
            assert!(res.is_start);
        }

        // case 2: +/+ match at end
        // flag 0 (+/+), 144M6S => the read aligned to the top strand of the transposon, and the match was at the end of the transposon
        //           match ... ---->
        //              -------->
        //              <--------
        // expected: m_size = 144, s_size = 6, is_sm = false, is_start = false
        let input2 = "2L_Read_977219	0	roo#LTR/Bel-Pao	8949	0	144M6S	*	0	0	GGACTATTTACGTAGGCCTCTGCGTAGGCCATTTACTTTAAGATGCGATTCTCATGTCACCTATTTAAACCGAAGATATTTCCAAATAAAACCAGTTTCTTACAAAAACTCAACGAGTAAAGTCTTCTTATTTGGGATTTTACATTTGGT	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	NM:i:2	MD:Z:91T36C15	AS:i:134	XS:i:133	XA:Z:roo#LTR/Bel-Pao,+285,95M1D55M,3;";
        let result2 = TeAlignment::create(String::from(input2), &transposons_arc);
        assert!(result2.is_some());
        if let Some(res) = result2 {
            assert_eq!(res.m_size, 144);
            assert_eq!(res.s_size, 6);
            assert!(!res.is_sm);
            assert!(!res.is_start);
        }

        // case 3: +/- match at start
        // flag 16 (+/-), 9S141M => the read aligned to the bottom strand of the transposon, and the match was at the beginning of the transposon:
        //              -------->
        //              <--------
        // match ... <----
        // expected: m_size = 141, s_size = 9, is_sm = true, is_start = true
        let input3 = "2L_Read_355243	16	blood#LTR/Gypsy	1	0	9S141M	*	0	0	GTGGCGAATTGTAGTATGTGCATATATCGAGGGTATACTGTACCTATAAGTACACAGCAACACTTAGTTGCATTGCATAAATAAATGTCTCAAGTGAGCGTGATATAAGATCACCCATTTATGCTTTAAGCTAAGTCAGCATCCCCACGC	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	NM:i:1	MD:Z:26C114	AS:i:136	XS:i:136	XA:Z:blood#LTR/Gypsy,-7012,9S141M,1;";
        let result3 = TeAlignment::create(String::from(input3), &transposons_arc);
        assert!(result3.is_some());
        if let Some(res) = result3 {
            assert_eq!(res.m_size, 141);
            assert_eq!(res.s_size, 9);
            assert!(res.is_sm);
            assert!(res.is_start);
        }

        // case 4: +/- match at end
        // flag 16 (+/-), 31M119S => the read aligned to the top strand of the transposon, and the match was at the end of the transposon
        //              -------->
        //              <--------
        //          match ... <----
        // expected: m_size = 31, s_size = 119, is_sm = false, is_start = false
        let input4 = "2L_Read_347822	16	blood#LTR/Gypsy	7380	0	31M119S	*	0	0	CTCAATTGGTGGCATATATTGGTTTATTACAGAATATCGAATCACTGATTCGGGATGTGAGAGTCACAATTTATTCCGCGATATCAGTTAAAAAAAATCTTCAAGACTTAAGATTTGACCGACAAAGAACATTTCTACGTGTTGGCCAAG	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	NM:i:0	MD:Z:31	AS:i:31	XS:i:31	XA:Z:blood#LTR/Gypsy,-368,31M119S,0;";
        let result4 = TeAlignment::create(String::from(input4), &transposons_arc);
        assert!(result4.is_some());
        if let Some(res) = result4 {
            assert_eq!(res.m_size, 31);
            assert_eq!(res.s_size, 119);
            assert!(!res.is_sm);
            assert!(!res.is_start);
        }
    }
}