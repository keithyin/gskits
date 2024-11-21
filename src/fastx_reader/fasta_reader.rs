use std::{
    fs::File,
    io::{self, BufRead, BufReader},
};

use super::{fastx_header_line_to_header, ReadInfo};

pub struct FastaFileReader {
    fname: String,
    lines: Option<io::Lines<BufReader<File>>>,
    current_header: Option<String>,
}

impl FastaFileReader {
    pub fn new(fname: String) -> Self {
        Self {
            fname: fname,
            lines: None,
            current_header: None,
        }
    }
    pub fn get_fname(&self) -> &str {
        return &self.fname;
    }
}

impl Iterator for FastaFileReader {
    type Item = ReadInfo;
    fn next(&mut self) -> Option<Self::Item> {
        if self.lines.is_none() {
            let file = File::open(&self.fname).expect(&format!("open file error: {}", self.fname));
            let lines = BufReader::new(file).lines();
            self.lines = Some(lines);
            let first_header = self.lines.as_mut().unwrap().next().unwrap().unwrap();
            assert!(first_header.starts_with(">"));

            self.current_header = Some(fastx_header_line_to_header(&first_header));
        }

        if self.current_header.is_none() {
            return None;
        }

        let mut fasta_record =
            ReadInfo::new_fa_record(self.current_header.take().unwrap(), String::new());

        while let Some(line) = self.lines.as_mut().unwrap().next() {
            let line = line.unwrap();
            if line.starts_with(">") {
                self.current_header = Some(fastx_header_line_to_header(&line));
                return Some(fasta_record);
            }

            fasta_record.seq.push_str(line.trim());
        }

        return Some(fasta_record);
    }
}
