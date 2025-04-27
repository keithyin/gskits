use std::collections::{HashMap, HashSet};

use crate::gsbam::bam_record_ext::{BamRecord, BamRecordExt};
pub mod region;
/// common data structures

#[derive(Debug, Default)]
pub struct ReadInfo {
    pub name: String,
    pub seq: String,
    pub cx: Option<u8>,
    pub ch: Option<u32>,
    pub np: Option<u32>,
    pub rq: Option<f32>,
    pub qual: Option<Vec<u8>>, // phreq, no offset
    pub dw: Option<Vec<u8>>,
    pub ar: Option<Vec<u8>>,
    pub cr: Option<Vec<u8>>,
    pub be: Option<Vec<u32>>,
    pub nn: Option<Vec<u8>>,
    pub wd: Option<Vec<u8>>, // linker width
    pub sd: Option<Vec<u8>>, // standard devition
    pub sp: Option<Vec<u8>>, // slope
}

impl ReadInfo {
    pub fn new_fa_record(name: String, seq: String) -> Self {
        let mut res = Self::default();
        res.name = name;
        res.seq = seq;
        res
    }

    pub fn new_fq_record(name: String, seq: String, qual: Vec<u8>) -> Self {
        let mut res = ReadInfo::new_fa_record(name, seq);
        res.qual = Some(qual);
        res
    }

    pub fn from_bam_record(
        record: &BamRecord,
        qname_suffix: Option<&str>,
        tags: &HashSet<String>,
    ) -> Self {
        let mut qname = unsafe { String::from_utf8_unchecked(record.qname().to_vec()) };
        if let Some(suffix) = qname_suffix {
            qname.push_str(suffix);
        }
        let record_ext = BamRecordExt::new(record);
        let seq = unsafe { String::from_utf8_unchecked(record.seq().as_bytes()) };

        let dw = if tags.contains("dw") {
            record_ext
                .get_dw()
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        let ar = if tags.contains("ar") {
            record_ext
                .get_ar()
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        let cr = if tags.contains("cr") {
            record_ext
                .get_cr()
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        let nn = if tags.contains("nn") {
            record_ext
                .get_nn()
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        let wd = if tags.contains("wd") {
            record_ext
                .get_uint_list(b"wd")
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        let sd = if tags.contains("sd") {
            record_ext
                .get_uint_list(b"sd")
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        let sp = if tags.contains("sp") {
            record_ext
                .get_uint_list(b"sd")
                .map(|v| v.into_iter().map(|v| v as u8).collect())
        } else {
            None
        };

        Self {
            name: qname,
            seq: seq,
            cx: record_ext.get_cx(),
            ch: record_ext.get_ch(),
            np: record_ext.get_np().map(|v| v as u32),
            rq: record_ext.get_rq(),
            qual: Some(record_ext.get_qual().to_vec()),
            dw: dw,
            ar: ar,
            cr: cr,
            be: record_ext.get_be(),
            nn: nn,
            wd: wd,
            sd: sd,
            sp: sp,
        }
    }
}

pub fn name2idx_and_seq(read_infos: &Vec<ReadInfo>) -> HashMap<&str, (usize, &str)> {
    read_infos
        .iter()
        .enumerate()
        .map(|(idx, read_info)| (read_info.name.as_str(), (idx, read_info.seq.as_str())))
        .collect()
}

pub fn name2idx(read_infos: &Vec<ReadInfo>) -> HashMap<&str, usize> {
    read_infos
        .iter()
        .enumerate()
        .map(|(idx, read_info)| (read_info.name.as_str(), idx))
        .collect()
}

pub fn name2seq(read_infos: &Vec<ReadInfo>) -> HashMap<&str, &str> {
    read_infos
        .iter()
        .map(|read_info| (read_info.name.as_str(), read_info.seq.as_str()))
        .collect()
}

pub fn idx2name_and_seq(read_infos: &Vec<ReadInfo>) -> HashMap<usize, (&str, &str)> {
    read_infos
        .iter()
        .enumerate()
        .map(|(idx, read_info)| (idx, (read_info.name.as_str(), read_info.seq.as_str())))
        .collect()
}

#[cfg(test)]
mod test {

    #[test]
    fn test_read_info() {}
}
