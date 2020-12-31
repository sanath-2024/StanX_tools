use path_abs::PathFile;

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::str::FromStr;

// represent insertion locations (within the reference)
// and read them in from a TSV file
#[derive(PartialEq, Eq, Hash)]
pub struct ILoc {
    pub chrom: String,
    pub upstream_pos: u64,
    pub downstream_pos: u64,
}

impl ILoc {
    pub fn length(&self) -> u64 {
        self.downstream_pos - self.upstream_pos + 1
    }
    fn read_line(line: String) -> Option<ILoc> {
        let fields: Vec<&str> = (&line[..]).split("\t").collect();
        if fields[7] == "reference" {
            return Some(ILoc {
                chrom: fields[0].to_owned(),
                upstream_pos: FromStr::from_str(fields[1]).unwrap(),
                downstream_pos: FromStr::from_str(fields[2]).unwrap(),
            });
        } else {
            return None;
        }
    }
    pub fn read_file(file_path: PathFile) -> Vec<ILoc> {
        let mut res = Vec::new();
        let lines = BufReader::new(File::open(file_path).unwrap())
            .lines()
            .map(|l| l.unwrap());
        let mut header_line = true;
        for line in lines {
            if header_line {
                header_line = false;
                continue;
            }
            if let Some(iloc) = ILoc::read_line(line) {
                // insertion locations should already be in the correct order:
                // sorted first by chromosome, then by upstream position
                res.push(iloc);
            }
        }
        return res;
    }
    pub fn contains(&self, chrom: &String, pos: u64) -> bool {
        pos >= self.upstream_pos && pos <= self.downstream_pos && chrom == &self.chrom
    }
}
