use lazy_static::lazy_static;
use regex::Regex;

// lazy static regexes to avoid multiple compilations
lazy_static! {
    pub static ref SM_REGEX: Regex = Regex::new(r"^(\d+)S(\d+)M$").unwrap();
    pub static ref MS_REGEX: Regex = Regex::new(r"^(\d+)M(\d+)S$").unwrap();
    pub static ref HM_REGEX: Regex = Regex::new(r"^(\d+)H(\d+)M$").unwrap();
    pub static ref MH_REGEX: Regex = Regex::new(r"^(\d+)M(\d+)H$").unwrap();
    pub static ref M_REGEX: Regex = Regex::new(r"^\d+M$").unwrap();
}

// get a numeric capture from a regex less verbosely
// steps:
// 1. unwrap the capture (assume that it exists)
// 2. get the Match object of the capture in position "which"
// 3. unwrap (assume that it exists, since the Regex has already been checked)
// 4. turn the Match into an &str
// 5. parse out the u64
// 6. unwrap (assume that it is a number, which we know because of the Regex itself)
pub fn get_capture<'t>(captures: Option<regex::Captures<'t>>, which: usize) -> u64 {
    captures
        .unwrap()
        .get(which)
        .unwrap()
        .as_str()
        .parse()
        .unwrap()
}
