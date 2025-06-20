use std::fmt::Display;

use rust_htslib::bam::{ext::BamRecordExtensions, record::Aux, record::Cigar, Record};

pub type BamRecord = rust_htslib::bam::Record;
pub type BamWriter = rust_htslib::bam::Writer;
pub type BamReader = rust_htslib::bam::Reader;

pub fn record2str(record: &Record) -> String {
    format!(
        "refstart:{}, refend:{}, seq:{}",
        record.reference_start(),
        record.reference_end(),
        String::from_utf8(record.seq().as_bytes()).unwrap()
    )
}

pub struct BamRecordExt<'a> {
    bam_record: &'a BamRecord,
    qname: Option<String>,
    seq: Option<String>,
}

impl<'a> BamRecordExt<'a> {
    pub fn new(bam_record: &'a BamRecord) -> Self {
        BamRecordExt {
            bam_record,
            qname: None,
            seq: None,
        }
    }

    /// eq / (eq + diff + ins + del)
    pub fn compute_identity(&self) -> f32 {
        let mut aligned_span = 0;
        let mut matched = 0;
        self.bam_record
            .cigar()
            .iter()
            .for_each(|cigar| match *cigar {
                Cigar::Equal(n) => {
                    matched += n;
                    aligned_span += n;
                }
                Cigar::Diff(n) => aligned_span += n,
                Cigar::Del(n) | Cigar::Ins(n) => aligned_span += n,
                Cigar::Match(_) => panic!("must use eqx cigar"),
                _ => {}
            });

        aligned_span = if aligned_span > 0 { aligned_span } else { 1 };

        matched as f32 / aligned_span as f32
    }

    /// (eq + diff + ins) / query_len
    pub fn compute_query_coverage(&self) -> f32 {
        let mut seq_len = self.bam_record.seq_len_from_cigar(true);
        let aligned_span = self
            .bam_record
            .cigar()
            .iter()
            .map(|cigar| match *cigar {
                Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Match(n) | Cigar::Ins(n) => n,
                _ => 0,
            })
            .reduce(|acc, b| acc + b)
            .unwrap_or(0);

        seq_len = if seq_len > 0 { seq_len } else { 1 };

        aligned_span as f32 / seq_len as f32
    }

    pub fn query_alignment_start(&self) -> usize {
        match self.bam_record.cigar().first() {
            Some(cigar) => match *cigar {
                Cigar::SoftClip(n) => n as usize,
                _ => 0,
            },

            _ => 0,
        }
    }

    pub fn query_alignment_end(&self) -> usize {
        let aligned_qlen = self
            .bam_record
            .cigar()
            .iter()
            .map(|cigar| match *cigar {
                Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) | Cigar::Match(n) => n,
                _ => 0,
            })
            .sum::<u32>();

        self.query_alignment_start() + aligned_qlen as usize
    }

    pub fn reference_start(&self) -> usize {
        assert!(self.bam_record.pos() >= 0, "set_pos first");
        self.bam_record.pos() as usize
    }

    pub fn reference_end(&self) -> usize {
        assert!(self.bam_record.pos() >= 0, "set_pos first");
        let aligned_rlen = self
            .bam_record
            .cigar()
            .iter()
            .map(|cigar| match *cigar {
                Cigar::Equal(n)
                | Cigar::Diff(n)
                | Cigar::Del(n)
                | Cigar::Match(n)
                | Cigar::RefSkip(n) => n,
                _ => 0,
            })
            .sum::<u32>();

        self.bam_record.pos() as usize + aligned_rlen as usize
    }

    pub fn get_seq_cached(&mut self) -> &str {
        if self.seq.is_none() {
            self.seq = Some(self.get_seq());
        }

        self.seq.as_ref().unwrap()
    }

    pub fn get_qname_cached(&mut self) -> &str {
        if self.qname.is_none() {
            self.qname = Some(self.get_qname());
        }
        self.qname.as_ref().unwrap()
    }

    pub fn get_seq(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.bam_record.seq().as_bytes()) }
    }

    pub fn get_qname(&self) -> String {
        unsafe { String::from_utf8_unchecked(self.bam_record.qname().to_vec()) }
    }

    pub fn get_qual(&self) -> &[u8] {
        self.bam_record.qual()
    }

    pub fn get_cx(&self) -> Option<u8> {
        self.get_int(b"cx").map(|v| v as u8)
    }

    pub fn get_ch(&self) -> Option<u32> {
        self.get_int(b"ch")
    }

    pub fn get_np(&self) -> Option<u32> {
        self.get_int(b"np")
    }

    pub fn get_coverage(&self) -> Option<f32> {
        self.get_float(b"ec")
    }

    pub fn get_identity(&self) -> Option<f32> {
        self.get_float(b"iy")
    }

    pub fn get_cq(&self) -> Option<f32> {
        self.get_float(b"cq")
    }

    pub fn get_oe(&self) -> Option<f32> {
        self.get_float(b"oe")
    }

    pub fn get_rq(&self) -> Option<f32> {
        self.get_float(b"rq")
    }

    pub fn get_dw(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"dw")
    }

    pub fn get_ar(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"ar")
    }

    pub fn get_cr(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"cr")
    }

    pub fn get_float_cr(&self) -> Option<Vec<f32>> {
        self.get_float_list(b"cr")
    }

    pub fn get_nn(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"nn")
    }

    pub fn get_be(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"be")
    }

    pub fn get_str(&self, tag: &[u8]) -> Option<String> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::String(v) => Some(v.to_string()),
            _ => None,
        })
    }

    pub fn get_int(&self, tag: &[u8]) -> Option<u32> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::I8(v) => Some(v as u32),
            Aux::U8(v) => Some(v as u32),
            Aux::I16(v) => Some(v as u32),
            Aux::U16(v) => Some(v as u32),
            Aux::I32(v) => Some(v as u32),
            Aux::U32(v) => Some(v as u32),

            _ => None,
        })
    }

    #[allow(unused)]
    pub fn get_uint(&self, tag: &[u8]) -> Option<u32> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::U8(v) => Some(v as u32),
            Aux::U16(v) => Some(v as u32),
            Aux::U32(v) => Some(v as u32),
            _ => None,
        })
    }

    pub fn get_float(&self, tag: &[u8]) -> Option<f32> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::Float(v) => Some(v as f32),
            Aux::Double(v) => Some(v as f32),
            _ => None,
        })
    }

    pub fn get_float_list(&self, tag: &[u8]) -> Option<Vec<f32>> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::ArrayFloat(v) => Some(v.iter().collect::<Vec<f32>>()),
            _ => None,
        })
    }

    pub fn get_uint_list(&self, tag: &[u8]) -> Option<Vec<u32>> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::ArrayU8(v) => Some(v.iter().map(|v| v as u32).collect::<Vec<u32>>()),
            Aux::ArrayU16(v) => Some(v.iter().map(|v| v as u32).collect::<Vec<u32>>()),
            Aux::ArrayU32(v) => Some(v.iter().map(|v| v as u32).collect::<Vec<u32>>()),
            _ => None,
        })
    }

    #[allow(unused)]
    pub fn get_int_list(&self, tag: &[u8]) -> Option<Vec<i32>> {
        self.bam_record.aux(tag).ok().and_then(|aux| match aux {
            Aux::ArrayI8(v) => Some(v.iter().map(|v| v as i32).collect::<Vec<i32>>()),
            Aux::ArrayI16(v) => Some(v.iter().map(|v| v as i32).collect::<Vec<i32>>()),
            Aux::ArrayI32(v) => Some(v.iter().map(|v| v as i32).collect::<Vec<i32>>()),
            _ => None,
        })
    }
}

impl<'a> Display for BamRecordExt<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "refstart:{}, refend:{}, seq:{}",
            self.bam_record.reference_start(),
            self.bam_record.reference_end(),
            self.get_seq()
        )
    }
}

pub fn draw_aligned_seq(
    record: &BamRecord,
    ref_seq: &[u8],
    r_start: Option<usize>,
    r_end: Option<usize>,
) -> (String, String) {
    let mut ref_aligned_seq = String::new();
    let mut query_aligned_seq = String::new();

    let query_seq = record.seq().as_bytes();
    let mut rpos_cursor = None;
    for [qpos, rpos] in record.aligned_pairs_full() {
        if rpos.is_some() {
            rpos_cursor = rpos;
        }

        if let Some(r_start) = r_start {
            if let Some(rpos_cursor) = rpos_cursor {
                if (rpos_cursor as usize) < r_start {
                    continue;
                }
            } else {
                continue;
            }
        }

        let q_char = if let Some(qpos_) = qpos {
            unsafe { (*query_seq.get_unchecked(qpos_ as usize)) as char }
        } else {
            '-'
        };

        let r_char = if let Some(rpos_) = rpos {
            unsafe { (*ref_seq.get_unchecked(rpos_ as usize)) as char }
        } else {
            '-'
        };

        ref_aligned_seq.push(r_char);
        query_aligned_seq.push(q_char);

        if let Some(r_end) = r_end {
            if let Some(rpos_cursor) = rpos_cursor {
                if (rpos_cursor as usize) >= (r_end - 1) {
                    break;
                }
            }
        }
    }

    (ref_aligned_seq, query_aligned_seq)
}
