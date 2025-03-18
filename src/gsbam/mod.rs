use std::collections::HashSet;

use bam_record_ext::BamReader;
use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView, Read};

use crate::ds::ReadInfo;

pub mod bam_reader;
pub mod bam_record_ext;
pub mod cigar_ext;
pub mod bam_header_ext;
pub mod plp_counts_from_records;
pub mod query_locus_blacklist_gen;
pub mod utils;

#[deprecated(since="0.10.0", note="use gsbam::bam_header_ext::BamHeaderExt instead")]
pub fn get_last_pg_from_bam_header(header_view: &HeaderView) -> Option<String> {
    let header = Header::from_template(header_view);
    let header = header.to_hashmap();
    if let Some(pg_info) = header.get("PG") {
        if let Some(last) = pg_info.last() {
            return Some(
                last.get("ID")
                    .expect(&format!("No ID in PG header"))
                    .to_string(),
            );
        } else {
            return None;
        }
    } else {
        return None;
    }
}

static PG: &'static str = "PG";
pub fn build_pg_header(
    id: &str,
    pn: &str,
    cl: &str,
    vn: &str,
    pp: Option<&str>,
) -> HeaderRecord<'static> {
    let mut hd = HeaderRecord::new(PG.as_bytes());
    hd.push_tag(b"ID", id)
        .push_tag(b"PN", pn)
        .push_tag(b"CL", cl)
        .push_tag(b"VN", vn);
    if let Some(pp_) = pp {
        hd.push_tag(b"PP", pp_);
    }
    hd
}

pub fn read_bam(bam_file: &str, threads: Option<usize>) -> Vec<ReadInfo> {
    let threads = threads.unwrap_or(4);
    let mut reader = BamReader::from_path(bam_file).unwrap();
    reader.set_threads(threads).unwrap();
    
    reader
        .records()
        .into_iter()
        .map(|record| record.unwrap())
        .map(|record| ReadInfo::from_bam_record(&record, None, &HashSet::new()))
        .collect()
}
