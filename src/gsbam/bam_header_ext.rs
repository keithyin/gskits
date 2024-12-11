use rust_htslib::bam::{Header, HeaderView};

#[derive(Debug, Clone)]
pub struct HeaderSQ {
    tid: i32,
    name: String,
    len: usize,
}

impl HeaderSQ {
    pub fn new(tid: i32, name: String, len: usize) -> Self {
        Self { tid, name, len }
    }

    pub fn get_tid(&self) -> i32 {
        self.tid
    }

    pub fn get_name(&self) -> &str {
        &self.name
    }

    pub fn get_len(&self) -> usize {
        self.len
    }
}

pub struct BamHeaderExt {
    header: Header,
    last_pg: Option<String>,
    all_seqs: Option<Vec<HeaderSQ>>,
}

impl BamHeaderExt {
    pub fn new(header: Header) -> Self {
        Self {
            header,
            last_pg: None,
            all_seqs: None,
        }
    }

    pub fn get_last_pg_from_bam_header_cached(&mut self) -> Option<&str> {
        if self.last_pg.is_none() {
            let header = self.header.to_hashmap();
            if let Some(pg_info) = header.get("PG") {
                if let Some(last) = pg_info.last() {
                    self.last_pg = Some(
                        last.get("ID")
                            .expect(&format!("No ID in PG header"))
                            .to_string(),
                    );
                }
            }
        }

        self.last_pg.as_ref().map(|v| v.as_str())
    }

    pub fn get_all_seqs_cached(&mut self) -> Option<&Vec<HeaderSQ>> {
        if self.all_seqs.is_none() {
            let header = self.header.to_hashmap();
            if let Some(seq_infos) = header.get("SQ") {
                let mut header_seqs = vec![];
                seq_infos.iter().enumerate().for_each(|(tid, seq_info)| {
                    let name = seq_info.get("SN").unwrap().clone();
                    let len = seq_info.get("LN").unwrap().parse::<usize>().unwrap();
                    header_seqs.push(HeaderSQ::new(tid as i32, name, len));
                });

                if header_seqs.len() > 0 {
                    self.all_seqs = Some(header_seqs);
                }
            }
        }

        self.all_seqs.as_ref()
    }
}

impl From<&HeaderView> for BamHeaderExt {
    fn from(value: &HeaderView) -> Self {
        BamHeaderExt::new(Header::from_template(value))
    }
}
