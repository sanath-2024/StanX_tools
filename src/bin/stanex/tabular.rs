use anyhow::{Context, Result};

use std::collections::HashMap;

pub struct Metadata {
    pub delimiter: String,
    // a map between positions and headings
    // note: the positions are 1-indexed to make it
    // easy to look at the file format and generate this struct
    // also, not all headings in the file need to be present in this map ...
    // only the important ones
    pub headings: HashMap<usize, String>,
}

pub struct Data {
    // key: heading name
    // value: the data for that heading
    fields: HashMap<String, String>,
}

impl Metadata {
    pub fn read(&self, row: String) -> Data {
        let split_str: Vec<&str> = row.split(&self.delimiter[..]).collect();
        let mut res = Data {
            fields: HashMap::new(),
        };
        for (position, heading) in &self.headings {
            if position > &split_str.len() {
                panic!(
                    "error reading tabular data: position {} is greater than the number of columns ({}) ... string: \"{}\"",
                    position, split_str.len(), row
                );
            }
            res.fields
                .insert(heading.clone(), split_str[position - 1].to_string());
        }
        return res;
    }
}

impl Data {
    pub fn get(&self, heading: &str) -> Result<String> {
        let value = self
            .fields
            .get(&heading.to_string())
            .context(format!("field {} is invalid", heading))?;
        Ok(value.clone())
    }
}
