use std::collections::{HashSet, HashMap};
use rust_htslib::bam::{self, Read};

use super::bam_record_ext::BamRecordExt;




#[allow(unused)]
pub fn seperate_dup_nondup_channels(bam_filepath: &str) -> (HashSet<usize>, HashSet<usize>) {
    let mut bam_file = bam::Reader::from_path(bam_filepath).unwrap();
    bam_file.set_threads(16).unwrap();
    let mut counter = HashMap::new();
    for record in bam_file.records() {
        let record = record.unwrap();
        let record_ext = BamRecordExt::new(&record);
        
        let ch = if let Some(n) = record_ext.get_ch() {
            n as usize
        } else {
            // 20240402_Sync_Y0003_02_H01_Run0002_called_subreads/85550/ccs  get ch from qname
            record_ext.get_qname().split("/")
                .collect::<Vec<&str>>()[1].parse::<usize>().unwrap()
        };

        * counter.entry(ch).or_insert(0) += 1;
    }

    let nondup_chs = counter.iter()
        .filter(|(_, v)| **v == 1)
        .map(|(k, _)| *k)
        .collect::<HashSet<usize>>();

    let dup_chs = counter.iter()
        .filter(|(_, v)| **v > 1)
        .map(|(k, _)| *k)
        .collect::<HashSet<usize>>();

    (dup_chs, nondup_chs)
}