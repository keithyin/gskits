use crate::gsbam::bam_record_ext::{BamRecord, BamRecordExt};

/// common data structures


pub struct ReadInfo {
    pub name: String,
    pub seq: String, 
    pub ch: Option<u32>,
    pub np: Option<u32>,
    pub rq: Option<f32>,
    pub qual: Option<Vec<u8>>, // phreq, no offset
    pub dw: Option<Vec<u8>>,
    pub ar: Option<Vec<u8>>,
    pub cr: Option<Vec<u8>>,
    pub be: Option<Vec<u32>>
    
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
            cr: None,
            be: None,
        }
    }

    pub fn new_fq_record(name: String, seq: String, qual: Vec<u8>) -> Self {
        Self {
            name: name,
            seq: seq,
            ch: None,
            np: None,
            rq: None,
            qual: Some(qual),
            dw: None,
            ar: None,
            cr: None,
            be: None,

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
            qual: Some(record_ext.get_qual().to_vec()),
            dw: record_ext.get_dw().map(|v| v.into_iter().map(|v| v as u8).collect()),
            ar: record_ext.get_ar().map(|v| v.into_iter().map(|v| v as u8).collect()),
            cr: record_ext.get_cr().map(|v| v.into_iter().map(|v| v as u8).collect()),
            be: record_ext.get_be(),
        }
    }

}