use core::fmt;
use std::{cmp, collections::HashMap};

use rust_htslib::bam::{ext::BamRecordExtensions, IndexedReader, Read};

use super::bam_record_ext::{BamRecord, BamRecordExt};

use lazy_static::lazy_static;

lazy_static! {
    static ref FWD_BASE_2_IDX: HashMap<u8, usize> = {
        let mut m = HashMap::new();
        m.insert('A' as u8, 4);
        m.insert('C' as u8, 5);
        m.insert('G' as u8, 6);
        m.insert('T' as u8, 7);
        m.insert(' ' as u8, 9);
        m.insert('-' as u8, 9);
        m.insert('*' as u8, 9);
        m
    };
    static ref REV_BASE_2_IDX: HashMap<u8, usize> = {
        let mut m = HashMap::new();
        m.insert('A' as u8, 0);
        m.insert('C' as u8, 1);
        m.insert('G' as u8, 2);
        m.insert('T' as u8, 3);
        m.insert(' ' as u8, 8);
        m.insert('-' as u8, 8);
        m.insert('*' as u8, 8);

        m
    };
}

// lazy_static! {
//     pub static ref FWD_BASE_2_IDX: HashMap<u8, usize> = {
//         let mut map = HashMap::new();
//         map.insert('A' as u8, 0);
//         map.insert('C' as u8, 1);
//         map.insert('G' as u8, 2);
//         map.insert('T' as u8, 3);
//         map.insert('-' as u8, 8);
//         map
//     };
//     pub static ref REV_BASE_2_IDX: HashMap<u8, usize> = {
//         let mut map = HashMap::new();
//         map.insert('A' as u8, 4);
//         map.insert('C' as u8, 5);
//         map.insert('G' as u8, 6);
//         map.insert('T' as u8, 7);
//         map.insert('-' as u8, 9);
//         map
//     };
// }

pub fn get_base_idx(base: u8, fwd: bool) -> usize {
    if fwd {
        *FWD_BASE_2_IDX.get(&base).unwrap()
    } else {
        *REV_BASE_2_IDX.get(&base).unwrap()
    }
}

#[derive(Clone)]
pub struct PlpCnts {
    ref_start: usize,
    ref_end: usize,
    major: Vec<usize>,
    minor: Vec<usize>,
    cnts: Vec<u32>, // 10 * lengths. atcgATCG gap GAP
    major_start_idx: HashMap<usize, usize>,
    timesteps: usize,
}

impl fmt::Debug for PlpCnts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PlpCnts")
            .field("ref_start", &self.ref_start)
            .field("ref_end", &self.ref_end)
            .field("major", &self.major)
            .field("minor", &self.minor)
            .field("cnts", &self.cnts)
            .finish()
    }
}

impl PlpCnts {
    /// ref_pos_length, the vec is sorted by ref_pos
    pub fn new(ref_pos_length: Vec<(usize, usize)>) -> Self {
        let tot_lengths = ref_pos_length.iter().map(|&(_, len)| len).sum::<usize>();
        let mut major = vec![0; tot_lengths];
        let mut minor = vec![0; tot_lengths];
        let mut major_start_idx = HashMap::new();
        let mut idx = 0;
        ref_pos_length.iter().for_each(|&(ref_pos, len)| {
            major_start_idx.insert(ref_pos, idx);
            for minor_pos in 0..len {
                assert!(idx < tot_lengths);
                unsafe {
                    *major.get_unchecked_mut(idx) = ref_pos;
                    *minor.get_unchecked_mut(idx) = minor_pos;
                }

                idx += 1;
            }
        });

        let ref_start = *major.first().unwrap();
        let ref_end = *major.last().unwrap() + 1;

        let timesteps = major.len();

        Self {
            ref_start,
            ref_end,
            major: major,
            minor: minor,
            cnts: vec![0; 10 * tot_lengths],
            major_start_idx,
            timesteps,
        }
    }

    pub fn update_with_records(&mut self, records: &Vec<BamRecord>) {
        for record in records {
            self.update_with_record(record);
        }
    }

    pub fn update_with_record(&mut self, record: &BamRecord) {
        let record_ext = BamRecordExt::new(record);
        let start = cmp::max(self.ref_start as i64, record_ext.reference_start() as i64);
        let end = cmp::min(self.ref_end as i64, record_ext.reference_end() as i64);
        let fwd = !record.is_reverse();
        let mut rpos_cursor = None;
        let mut qpos_cursor = None;
        let mut cur_ins = 0;
        let mut anchor = 0;
        let query_seq = record.seq().as_bytes();
        let query_end = record_ext.query_alignment_end();

        // println!("qname:{}, start:{}, end:{}", record_ext.get_qname(), start, end);
        for [qpos, rpos] in record.aligned_pairs_full() {
            if rpos.is_some() {
                rpos_cursor = rpos;
            }
            if qpos.is_some() {
                qpos_cursor = qpos;
            }
            if rpos_cursor.is_none() {
                continue;
            }

            if rpos_cursor.unwrap() < start {
                continue;
            }
            if rpos_cursor.unwrap() >= end {
                break;
            }

            if let Some(qpos_cursor_) = qpos_cursor {
                if qpos_cursor_ as usize >= query_end {
                    break;
                }
            }
            // print!("{},", rpos_cursor.unwrap());
            if let Some(rpos_) = rpos {
                if !self.major_start_idx.contains_key(&(rpos_ as usize)) {
                    let mut all_pos = self.major_start_idx.keys().map(|&v| v).collect::<Vec<_>>();
                    all_pos.sort();
                    panic!("rpos not found: {}, all pos are: {:?}", rpos_, all_pos);
                }

                anchor = *self.major_start_idx.get(&(rpos_ as usize)).unwrap();

                cur_ins = 0;
            } else {
                cur_ins += 1;
            }

            if let Some(qpos_) = qpos {
                self.update_cnts(anchor + cur_ins, query_seq[qpos_ as usize], fwd);
            } else {
                self.update_cnts(anchor + cur_ins, '-' as u8, fwd);
            }
        }
        // println!("");
    }

    pub fn get_major(&self) -> &Vec<usize> {
        &self.major
    }

    pub fn get_minor(&self) -> &Vec<usize> {
        &self.minor
    }

    pub fn get_cnts(&self) -> &Vec<u32> {
        &self.cnts
    }

    fn update_cnts(&mut self, col: usize, base: u8, fwd: bool) {
        let base_idx = get_base_idx(base, fwd);
        self.cnts[col + base_idx * self.timesteps] += 1;
    }

    pub fn cnts2str(&self) -> String {
        let mut tot_str_list = vec![];
        for col in 0..self.timesteps {
            let mut col_strs = vec![];
            for row in 0..10 {
                let idx = col + row * self.timesteps;
                col_strs.push(self.cnts[idx].to_string());
            }
            tot_str_list.push(col_strs.join(", "));
        }
        tot_str_list.join("\n")
    }
}

pub fn plp_within_region(
    reader: &mut IndexedReader,
    contig: &str,
    start: Option<usize>,
    end: Option<usize>,
) -> PlpCnts {
    if start.is_none() && end.is_none() {
        reader.fetch(contig).unwrap();
    }
    if start.is_none() ^ end.is_none() {
        panic!("start and end need to be all presented or all missed");
    }

    if start.is_some() {
        reader
            .fetch((contig, start.unwrap() as u32, end.unwrap() as u32))
            .unwrap();
    }

    let mut records = vec![];
    loop {
        let mut record = BamRecord::new();
        match reader.read(&mut record) {
            Some(res) => match res {
                Ok(_) => records.push(record),
                Err(_) => panic!("read record error"),
            },
            None => break,
        }
    }

    plp_with_records_region(&records, start, end)
}

pub fn plp_with_records_region(
    records: &Vec<BamRecord>,
    start: Option<usize>,
    end: Option<usize>,
) -> PlpCnts {
    let max_ins_of_ref_position = compute_max_ins_of_each_position(records, start, end);
    let len_of_ref_position = max_ins_of_ref_position
        .iter()
        .map(|(&pos, &ins)| (pos, ins + 1))
        .collect::<HashMap<_, _>>();

    let mut len_of_ref_positions_list = len_of_ref_position
        .iter()
        .map(|(&pos, &len)| (pos as usize, len as usize))
        .collect::<Vec<(_, _)>>();
    len_of_ref_positions_list.sort_by_key(|v| v.0);

    let mut plp_cnts = PlpCnts::new(len_of_ref_positions_list);
    plp_cnts.update_with_records(&records);
    plp_cnts
}

/// todo: use the cigar_str to speed up this function
pub fn compute_max_ins_of_each_position(
    records: &Vec<BamRecord>,
    rstart: Option<usize>,
    rend: Option<usize>,
) -> HashMap<i64, i32> {
    let mut pos2ins = HashMap::new();

    let rstart = rstart.map(|v| v as i64);
    let rend = rend.map(|v| v as i64);
    for record in records {
        let record_ext = BamRecordExt::new(record);
        let mut start = rstart.unwrap_or(record_ext.reference_start() as i64);
        let mut end = rend.unwrap_or(record_ext.reference_end() as i64);
        start = cmp::max(start, record_ext.reference_start() as i64);
        end = cmp::min(end, record_ext.reference_end() as i64);

        let mut rpos_cursor = None;
        // let mut qpos_cursor = None;
        let mut cur_ins = 0;
        let query_end = record_ext.query_alignment_end();

        // println!(
        //     "max_ins, qname:{}, rstart:{}, rend:{}",
        //     record_ext.get_qname(),
        //     start,
        //     end
        // );
        let mut aligned_pair_full = record.aligned_pairs_full().collect::<Vec<_>>();
        if aligned_pair_full.len() == 0 {
            continue;
        }

        // this make the following for loop correct.
        // if the query match to the last base of the ref seqence. the following for loop won't give the right result
        // but add this , the result is right
        if let Some(last_ref_pos) = aligned_pair_full.last().unwrap()[0] {
            aligned_pair_full.push([Some(last_ref_pos + 1), None]);
        }

        for [qpos, rpos] in aligned_pair_full.into_iter() {
            if rpos.is_some() {
                rpos_cursor = rpos;
            }
            // if qpos.is_some() {
            //     qpos_cursor = qpos;
            // }
            if rpos_cursor.is_none() {
                continue;
            }

            if rpos_cursor.unwrap() < start {
                continue;
            }

            // print!("{},", rpos_cursor.unwrap());
            if rpos_cursor.unwrap() >= end {
                // set the last rpos max ins and then break
                let rpos_ = rpos_cursor.unwrap();
                pos2ins.entry(rpos_ - 1).or_insert(0);
                *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                    cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                break;
            }

            if let Some(qpos_) = qpos {
                if qpos_ as usize >= query_end {
                    // println!("query hit end: {}", qpos_);
                    // set the last rpos max ins and then break
                    let rpos_ = rpos_cursor.unwrap();

                    pos2ins.entry(rpos_).or_insert(0);
                    *pos2ins.get_mut(&rpos_).unwrap() =
                        cmp::max(*pos2ins.get(&rpos_).unwrap(), cur_ins);
                    break;
                }
            }

            if let Some(rpos_) = rpos {
                if rpos_ > start {
                    pos2ins.entry(rpos_ - 1).or_insert(0);
                    *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                        cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                }

                cur_ins = 0;
            } else {
                cur_ins += 1;
            }
        }
        // println!("");
    }

    pos2ins
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use rust_htslib::bam::{ext::BamRecordExtensions, Header, IndexedReader, Read};

    use crate::gsbam::{
        bam_record_ext::{BamRecord, BamRecordExt},
        cigar_ext::parse_cigar_string,
        plp_counts_from_records::PlpCnts,
    };

    use super::{compute_max_ins_of_each_position, plp_within_region};

    #[test]
    fn test_test_plp_using_aligned_pairs_with_right_soft_clip() {
        // ACTC---
        // ACTCGGG

        let mut record = BamRecord::new();
        let seq = "ACTCGGG";
        println!("{}, {}", record.reference_start(), record.reference_end());
        println!("{}", i64::MAX);
        record.set_pos(0);

        record.set(
            b"qname",
            Some(&parse_cigar_string("4=3S").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );

        // NOTE: record.reference_end() does not work at this point. so i implemented reference_end in BamRecordExt

        let record_ext = BamRecordExt::new(&record);
        assert_eq!(record_ext.query_alignment_start(), 0);
        assert_eq!(record_ext.query_alignment_end(), 4);
        assert_eq!(record_ext.reference_start(), 0);
        assert_eq!(record_ext.reference_end(), 4);

        let position_max_ins = compute_max_ins_of_each_position(&vec![record], None, None);
        assert_eq!(position_max_ins.len(), 4);
        assert_eq!(position_max_ins.get(&0), Some(&0));
        assert_eq!(position_max_ins.get(&1), Some(&0));
        assert_eq!(position_max_ins.get(&2), Some(&0));
        assert_eq!(position_max_ins.get(&3), Some(&0));
    }

    #[test]
    fn test_test_plp_using_aligned_pairs_with_indel() {
        // --A-CTCC---
        // GGACCT-CGGG

        let mut record = BamRecord::new();
        let seq = "GGACCTCGGG";
        record.set_pos(0);

        record.set(
            b"qname",
            Some(&parse_cigar_string("2S1=1I2=1D1=3S").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );

        // NOTE: record.reference_end() does not work at this point. so i implemented reference_end in BamRecordExt

        let record_ext = BamRecordExt::new(&record);
        assert_eq!(record_ext.query_alignment_start(), 2);
        assert_eq!(record_ext.query_alignment_end(), 7);
        assert_eq!(record_ext.reference_start(), 0);
        assert_eq!(record_ext.reference_end(), 5);

        let position_max_ins = compute_max_ins_of_each_position(&vec![record], None, None);
        assert_eq!(position_max_ins.len(), 5);
        assert_eq!(position_max_ins.get(&0), Some(&1));
        assert_eq!(position_max_ins.get(&1), Some(&0));
        assert_eq!(position_max_ins.get(&2), Some(&0));
        assert_eq!(position_max_ins.get(&3), Some(&0));
        assert_eq!(position_max_ins.get(&4), Some(&0));
    }

    #[test]
    fn test_plp_cnts() {
        // --A--CTCC---
        // GGACCCT-CGGG
        //   A--TTCC
        //      CTCC

        let mut records = vec![];
        let mut record = BamRecord::new();
        let seq = "GGACCCTCGGG";
        record.set_pos(0);

        record.set(
            b"qname",
            Some(&parse_cigar_string("2S1=2I2=1D1=3S").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        records.push(record);

        let mut record = BamRecord::new();
        let seq = "ATTCC";
        record.set_pos(0);

        record.set(
            b"qname1",
            Some(&parse_cigar_string("1=1X3=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        records.push(record);

        let mut record = BamRecord::new();
        let seq = "CTCC";
        record.set_pos(1);

        record.set(
            b"qname2",
            Some(&parse_cigar_string("4=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        records.push(record);

        let max_ins_of_ref_position = compute_max_ins_of_each_position(&records, None, None);
        let len_of_ref_position = max_ins_of_ref_position
            .iter()
            .map(|(&pos, &ins)| (pos, ins + 1))
            .collect::<HashMap<_, _>>();

        let mut len_of_ref_positions_list = len_of_ref_position
            .iter()
            .map(|(&pos, &len)| (pos as usize, len as usize))
            .collect::<Vec<(_, _)>>();
        len_of_ref_positions_list.sort_by_key(|v| v.0);

        let mut plp_cnts = PlpCnts::new(len_of_ref_positions_list);
        plp_cnts.update_with_records(&records);
        assert_eq!(
            plp_cnts.get_cnts(),
            &vec![
                2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        );

        // println!("{}", plp_cnts.cnts2str());
    }

    #[test]
    fn test_region_plp_cnts() {
        // --A--CTCC---
        // GGACCCT-CGGG
        //   A--TTCC
        //      CTCC

        let mut records = vec![];
        let mut record = BamRecord::new();
        let seq = "GGACCCTCGGG";
        record.set_pos(0);

        record.set(
            b"qname0",
            Some(&parse_cigar_string("2S1=2I2=1D1=3S").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        records.push(record);

        let mut record = BamRecord::new();
        let seq = "ATTCC";
        record.set_pos(0);

        record.set(
            b"qname1",
            Some(&parse_cigar_string("1=1X3=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        records.push(record);

        let mut record = BamRecord::new();
        let seq = "CTCC";
        record.set_pos(1);

        record.set(
            b"qname2",
            Some(&parse_cigar_string("4=").unwrap()),
            seq.as_bytes(),
            &vec![255; seq.len()],
        );
        records.push(record);

        let max_ins_of_ref_position = compute_max_ins_of_each_position(&records, Some(1), Some(4));
        let len_of_ref_position = max_ins_of_ref_position
            .iter()
            .map(|(&pos, &ins)| (pos, ins + 1))
            .collect::<HashMap<_, _>>();

        let mut len_of_ref_positions_list = len_of_ref_position
            .iter()
            .map(|(&pos, &len)| (pos as usize, len as usize))
            .collect::<Vec<(_, _)>>();
        len_of_ref_positions_list.sort_by_key(|v| v.0);

        let mut plp_cnts = PlpCnts::new(len_of_ref_positions_list);
        plp_cnts.update_with_records(&records);
        // println!("{}", plp_cnts.cnts2str());
        // println!("{:?}", plp_cnts);
        assert_eq!(
            plp_cnts.get_cnts(),
            &vec![
                0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 0
            ]
        );
    }

    #[test]
    fn test_plp_within_region() {
        let mut reader = IndexedReader::from_path("test_data/sbr2smc.aligned.bam").unwrap();
        reader.set_threads(10).unwrap();

        let header = Header::from_template(reader.header());
        let header_hm = header.to_hashmap();
        let seqs = header_hm.get("SQ").unwrap();
        for seq in seqs {
            println!("{}", seq.get("SN").unwrap());
        }
        // plp_within_region(
        //     &mut reader,
        //     "20241108_Sync_Y0701_04_H01_Run0001/521756/smc",
        //     None,
        //     None,
        // );
    }
}
