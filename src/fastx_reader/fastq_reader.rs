use std::{
    fs::{self, File},
    io::BufReader,
};

use super::ReadInfo;
use bio::io::fastq::{self, FastqRead};

pub struct FastqReader {
    fname: String,
    reader: fastq::Reader<BufReader<File>>,
    record: fastq::Record,
}

impl FastqReader {
    pub fn new(fname: String) -> Self {
        let file = fs::File::open(&fname).expect(&format!("open {} error", fname));
        let reader = BufReader::new(file);
        let reader = fastq::Reader::from_bufread(reader);
        Self {
            fname,
            reader,
            record: fastq::Record::new(),
        }
    }
    pub fn get_fname(&self) -> &str {
        return &self.fname;
    }
}

impl Iterator for FastqReader {
    type Item = ReadInfo;
    fn next(&mut self) -> Option<Self::Item> {
        self.reader
            .read(&mut self.record)
            .expect(&format!("read fastq record error. {}", self.fname));

        if self.record.is_empty() {
            return None;
        }
        Some(ReadInfo::new_fq_record(
            self.record.id().to_string(),
            unsafe { String::from_utf8_unchecked(self.record.seq().to_vec()) },
            self.record.qual().iter().map(|v| *v - 33).collect(),
        ))

        // let qual = self
        //     .read_one_line()
        //     .expect("not a valid FastqRecord")
        //     .trim()
        //     .as_bytes()
        //     .iter()
        //     .map(|v| *v - 33)
        //     .collect()
        //     ;
    }
}
