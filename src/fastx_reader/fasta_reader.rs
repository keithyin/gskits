use std::{fs::File, io::BufReader};

use bio::io::fasta::{self, FastaRead};

use super::ReadInfo;

pub struct FastaFileReader {
    fname: String,
    reader: fasta::Reader<BufReader<File>>,
    record: fasta::Record,
}

impl FastaFileReader {
    pub fn new(fname: String) -> Self {
        let file = File::open(&fname).expect(&format!("open file error: {}", fname));
        let reader = BufReader::new(file);
        let reader = fasta::Reader::from_bufread(reader);

        Self {
            fname: fname,
            reader: reader,
            record: fasta::Record::new(),
        }
    }
    pub fn get_fname(&self) -> &str {
        return &self.fname;
    }
}

impl Iterator for FastaFileReader {
    type Item = ReadInfo;
    fn next(&mut self) -> Option<Self::Item> {
        self.reader
            .read(&mut self.record)
            .expect(&format!("read fasta error {}", self.fname));

        return if self.record.is_empty() {
            None
        } else {
            Some(ReadInfo::new_fa_record(
                self.record.id().to_string(),
                unsafe { String::from_utf8_unchecked(self.record.seq().to_vec()) },
            ))
        };
    }
}
