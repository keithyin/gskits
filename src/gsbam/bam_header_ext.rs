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

#[derive(Debug)]
pub struct PgHeader {
    pub id: Option<String>,
    pub program_name: Option<String>,
    pub program_version: Option<String>,
    pub prev_program: Option<String>,
    cmd_line: Option<String>,
}

impl PgHeader {
    pub fn get_cmd_line(&self) -> Option<String> {
        self.cmd_line
            .as_ref()
            .map(|line| line.split_whitespace().collect::<Vec<_>>().join(" "))
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

    pub fn get_last_pg_header(&self) -> Option<PgHeader> {
        let header = self.header.to_hashmap();
        return if let Some(pg_info) = header.get("PG") {
            if let Some(last) = pg_info.last() {
                Some(PgHeader {
                    id: last.get("ID").map(|item| item.to_string()),
                    program_name: last.get("PN").map(|item| item.to_string()),
                    program_version: last.get("VN").map(|item| item.to_string()),
                    prev_program: last.get("PP").map(|item| item.to_string()),
                    cmd_line: last.get("CL").map(|item| item.to_string()),
                })
            } else {
                None
            }
        } else {
            None
        };
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

#[cfg(test)]
mod test {
    use rust_htslib::bam::{Header, Read};

    use super::BamHeaderExt;

    #[test]
    fn test_header_ext() {
        let bam_file = rust_htslib::bam::Reader::from_path("test_data/header_only.bam").unwrap();
        let header = Header::from_template(bam_file.header());
        let header_ext = BamHeaderExt::new(header);
        println!("{:?}", header_ext.get_last_pg_header().map(|v| v.get_cmd_line()));
    }
}
