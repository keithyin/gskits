use std::{
    fs::{self, File},
    io::{BufRead, BufReader},
};

use super::{fastx_header_line_to_header, ReadInfo};

pub struct FastqReader {
    fname: String,
    reader: BufReader<File>
}

impl FastqReader {
    pub fn new(fname: String) -> Self {
        let file = fs::File::open(&fname).expect(&format!("open {} error", fname));
        let reader = BufReader::new(file);
        Self { fname, reader }
    }
    pub fn get_fname(&self) -> &str{
        return &self.fname;
    }

    fn read_one_line(&mut self) -> Option<String> {
        let mut line = String::new();
        if let Ok(n) = self.reader.read_line(&mut line) {
            if n == 0 {
                return None;
            }
            line = line.trim().to_string();
        } else {
            return None;
        }
        return Some(line);
    }
}

impl Iterator for FastqReader {
    type Item = ReadInfo;
    fn next(&mut self) -> Option<Self::Item> {
        let header = self.read_one_line();
        if header.is_none() || header.as_ref().unwrap().trim().len() == 0 {
            return None;
        }
        let mut header = header.unwrap();
        if !header.starts_with("@") {
            panic!("header:'{}' not a valid fastq header", header);
        }

        header = fastx_header_line_to_header(&header);

        let seq = self
            .read_one_line()
            .expect("not a valid FastqRecord")
            .trim()
            .to_string();
        let plus = self
            .read_one_line()
            .expect("not a valid FastqRecord")
            .trim()
            .to_string();

        if !plus.starts_with("+") {
            panic!("plus:'{}' not a valid fastq plus", plus);
        }

        let qual = self
            .read_one_line()
            .expect("not a valid FastqRecord")
            .trim()
            .to_string();

        Some(ReadInfo::new_fq_record(header, seq, qual))
    }
}
