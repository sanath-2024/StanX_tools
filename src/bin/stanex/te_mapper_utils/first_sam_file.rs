use anyhow::Result;
use lazy_static::lazy_static;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

use super::te_alignment::TeAlignment;
use crate::tabular::Metadata;

lazy_static! {
    static ref FIRST_SAM_FILE_TE_METADATA: Metadata = {
        let mut headings = HashMap::new();
        headings.insert(2, "TE_NAME".to_string());
        headings.insert(3, "TE_LEN".to_string());
        Metadata {
            delimiter: "\t".to_string(),
            headings: headings,
        }
    };
    static ref FIRST_SAM_FILE_ALIGNMENT_METADATA: Metadata = {
        let mut headings = HashMap::new();
        headings.insert(1, "QNAME".to_string());
        headings.insert(2, "FLAG".to_string());
        headings.insert(3, "RNAME".to_string());
        headings.insert(4, "POS".to_string());
        headings.insert(6, "CIGAR".to_string());
        headings.insert(10, "SEQ".to_string());
        Metadata {
            delimiter: "\t".to_string(),
            headings: headings,
        }
    };
}

fn read_te_into_map(te_str: String, transposon_lengths: &mut HashMap<String, u64>) {
    let te_data = FIRST_SAM_FILE_TE_METADATA.read(te_str);
    transposon_lengths.insert(
        te_data.get("TE_NAME").unwrap()[3..].to_string(),
        te_data.get("TE_LEN").unwrap()[3..].parse().unwrap(),
    );
}

pub fn read_all_tes_into_map(reader: &mut BufReader<File>) -> HashMap<String, u64> {
    // reads all TE's into a map and positions the buffered reader on the first line that is an alignment
    // read the file line by line
    // get rid of comments (comments in the SAM file start with "@SQ")
    // and ignore the last comment line (starts with "@PG")
    let mut transposon_lengths: HashMap<String, u64> = HashMap::new();

    let mut te_aligned_read;

    loop {
        te_aligned_read = String::new();
        reader.read_line(&mut te_aligned_read).unwrap();
        // get rid of trailing newline
        te_aligned_read = te_aligned_read[..te_aligned_read.len() - 1].to_string();
        if te_aligned_read.chars().nth(1).unwrap() == 'P' {
            break;
        } else {
            read_te_into_map(te_aligned_read, &mut transposon_lengths);
        }
    }

    return transposon_lengths;
}

pub fn read_te_alignment(
    alignment_str: String,
    transposon_lengths: &HashMap<String, u64>,
) -> Result<TeAlignment> {
    let alignment_data = FIRST_SAM_FILE_ALIGNMENT_METADATA.read(alignment_str);
    return TeAlignment::create(alignment_data, transposon_lengths);
}

#[cfg(test)]
mod tests {
    use super::{read_all_tes_into_map, read_te_alignment};

    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn test_cigar_string_parsing() {
        // test the creation of, m_size, s_size, is_sm, and is_start

        // read in the transposons
        let mut te_aligned_reader =
            BufReader::with_capacity(65_536, File::open("test/te_aligned.sam").unwrap());
        let transposon_lengths: HashMap<String, u64> =
            read_all_tes_into_map(&mut te_aligned_reader);

        // use actual inputs from the file
        // 54S34M62S is invalid
        let input0 = "I_MADE_UP_THIS_READ	0	roo#LTR/Bel-Pao	1	60	54S34M62S	*	0	0	CCTGGCTTGGGGCGGCCGCGGGTTCGTGGCGTCGGCGCTATTTGTTCCTTGGCAGTCGGCTCTTCCTATCATTGTGAAGCAAAATTCATATGGCATTGTCTCCTAAAACTTTTCTATAGTGCCGTATTTCTATGGCGCCCACTGTGAAGN	--F-7-F7----A--F7---------7--77----J7<---77--7---A--7-7-----<A7--7F<7FAAJA7---F-<-F7<<-<----<<--<---<--F7-F-F-<JFJAF7<JFJJAJFJFFJAJJJJJJJJJJJJFJJF<AA#	NM:i:0	MD:Z:34	AS:i:34	XS:i:0";
        let mut te_alignment = read_te_alignment(input0.to_string(), &transposon_lengths);
        assert!(te_alignment.is_err());

        // case 1: +/+ match at start
        // flag 0 (+/+), 119S31M => the read aligned to the top strand of the transposon, and the match was at the beginning of the transposon:
        //  match ... ---->
        //              -------->
        //              <--------
        // expected: m_size = 31, s_size = 119, is_sm = true, is_start = true
        let input1 = "2L_Read_976816	0	roo#LTR/Bel-Pao	1	0	119S31M	*	0	0	ACATATGATATAAATAGCATTAAATGTTGAGTATAACGTGTCAAAGAATCCTTGGGATGAATAATAACGGAGGAAGCTGTAAATATAACCAGATTAGAAACCTATTCCTATAAACTCTCTGTTCACACATGAACACGAATATATTTAAAG	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	NM:i:0	MD:Z:31	AS:i:31	XS:i:31	XA:Z:roo#LTR/Bel-Pao,+8665,119S31M,0;";
        te_alignment = read_te_alignment(input1.to_string(), &transposon_lengths);
        assert!(te_alignment.is_ok());
        if let Ok(res) = te_alignment {
            assert_eq!(res.qname, "2L_Read_976816");
            assert_eq!(res.rname, "roo#LTR/Bel-Pao");
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
        te_alignment = read_te_alignment(input2.to_string(), &transposon_lengths);
        assert!(te_alignment.is_ok());
        if let Ok(res) = te_alignment {
            assert_eq!(res.qname, "2L_Read_977219");
            assert_eq!(res.rname, "roo#LTR/Bel-Pao");
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
        te_alignment = read_te_alignment(input3.to_string(), &transposon_lengths);
        assert!(te_alignment.is_ok());
        if let Ok(res) = te_alignment {
            assert_eq!(res.qname, "2L_Read_355243");
            assert_eq!(res.rname, "blood#LTR/Gypsy");
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
        te_alignment = read_te_alignment(input4.to_string(), &transposon_lengths);
        assert!(te_alignment.is_ok());
        if let Ok(res) = te_alignment {
            assert_eq!(res.qname, "2L_Read_347822");
            assert_eq!(res.rname, "blood#LTR/Gypsy");
            assert_eq!(res.m_size, 31);
            assert_eq!(res.s_size, 119);
            assert!(!res.is_sm);
            assert!(!res.is_start);
        }
    }
}
