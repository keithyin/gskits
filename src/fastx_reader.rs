use crate::gsbam::bam_record_ext::{BamRecord, BamRecordExt};

pub mod fasta_reader;
pub mod fastq_reader;
pub mod fastx2bam;

/// 表示一个 FASTA 记录的结构体
#[derive(Debug, PartialEq)]
pub struct ReadsInfo {
    pub name: String,
    pub sequence: String,
    pub qual: Option<String>,  // phred + 33 string
    pub ch: Option<usize>,
    pub np: Option<usize>,
}

impl ReadsInfo {
    pub fn new_fa_record(name: String, seq: String) -> Self {
        Self {
            name: name,
            sequence: seq,
            qual: None,
            ch: None,
            np: None,
        }
    }

    pub fn new_fq_record(name: String, seq: String, qual: String) -> Self {
        Self {
            name: name,
            sequence: seq,
            qual: Some(qual),
            ch: None,
            np: None,
        }
    }

    pub fn from_bam_record(record: &BamRecord, qname_suffix: Option<&str>) -> Self {
        let mut qname = unsafe { String::from_utf8_unchecked(record.qname().to_vec()) };
        if let Some(suffix) = qname_suffix {
            qname.push_str(suffix);
        }

        let record_ext = BamRecordExt::new(record);

        let seq = unsafe { String::from_utf8_unchecked(record.seq().as_bytes()) };

        Self {
            name: qname,
            sequence: seq,
            qual: None,
            ch: record_ext.get_ch(),
            np: record_ext.get_np(),
        }
    }
}

pub fn fastx_header_line_to_header(header_line: &str) -> String {
    let header = if header_line[1..].contains(" ") {
        header_line[1..].split_once(" ").unwrap().0.to_string()
    } else {
        header_line[1..].to_string()
    };
    return header;
}


pub fn read_fastx<T: Iterator<Item = ReadsInfo>>(reader: T) -> Vec<ReadsInfo> {
    reader.into_iter().collect::<_>()
}