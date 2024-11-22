use crate::gsbam::bam_record_ext::{BamRecord, BamRecordExt};

/// common data structures


pub struct ReadInfo {
    pub name: String,
    pub seq: String, 
    pub ch: Option<usize>,
    pub np: Option<u32>,
    pub rq: Option<f32>,
    pub qual: Option<String>, // phreq+33 string
    pub dw: Option<u8>,
    pub ar: Option<u8>,
    pub cr: Option<u8>
}


impl ReadInfo {
    pub fn new_fa_record(name: String, seq: String) -> Self {
        Self {
            name: name,
            seq: seq,
            ch: None,
            np: None,
            rq: None,
            qual: None,
            dw: None,
            ar: None,
            cr: None
        }
    }

    pub fn new_fq_record(name: String, seq: String, qual: String) -> Self {
        Self {
            name: name,
            seq: seq,
            ch: None,
            np: None,
            rq: None,
            qual: Some(qual),
            dw: None,
            ar: None,
            cr: None

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
            seq: seq,
            ch: record_ext.get_ch(),
            np: record_ext.get_np().map(|v| v as u32),
            rq: record_ext.get_rq(),
            qual: None,
            dw: None,
            ar: None,
            cr: None

        }
    }
}