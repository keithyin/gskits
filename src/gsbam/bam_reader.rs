// use rust_htslib::bam::{Read, Records};

// use crate::fastx_reader::ReadsInfo;

// use super::bam_record_ext::BamReader;


// it's a self reference struct, hard to implement it using safe rust
// pub struct BamRecordReader<'a> {
//     fname: String,
//     bam_h: BamReader,
//     records_iter: Option<Records<'a, BamReader>>
// }

// impl <'a> BamRecordReader<'a> {
//     pub fn new(fname: String, threads: Option<usize>) -> Self {
//         let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
//         let mut bam_h = BamReader::from_path(&fname).expect(&format!("open {} error", fname));
//         bam_h.set_threads(threads).unwrap();
//         let mut out = BamRecordReader{fname, bam_h, records_iter: None};
//         out.records_iter = Some(out.bam_h.records());
//         out

//     }
// }

// impl<'a> Iterator for BamRecordReader<'a> {
//     type Item = ReadsInfo;
//     fn next(&mut self) -> Option<Self::Item> {
        
//         let mut bam_h = BamReader::from_path(&self.fname).expect(&format!("open {} error", self.fname));
//         let records: rust_htslib::bam::Records<'_, rust_htslib::bam::Reader> = bam_h.records();


//     }
// }