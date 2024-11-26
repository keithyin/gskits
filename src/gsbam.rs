use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView};

pub mod bam_reader;
pub mod bam_record_ext;

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
