use rust_htslib::bam::ext::BamRecordExtensions;

use crate::gsbam::bam_record_ext::{BamRecord, BamRecordExt};

#[derive(Debug, PartialEq, Eq)]
pub struct RefPolyLocusInfo {
    pub rstart: usize,
    pub rend: usize,
    pub qstart: usize,
    pub qend: usize,
    pub ref_base: char,
    pub ref_repeats: usize,
    pub query_repeats: usize,
    pub qseq: String,
    pub query_clean: bool,
}

impl RefPolyLocusInfo {
    fn new(
        rstart: usize,
        rend: usize,
        qstart: usize,
        qend: usize,
        ref_base: char,
        ref_repeats: usize,
        qseq: String,
    ) -> Self {
        let query_repeats = qseq
            .as_bytes()
            .iter()
            .filter(|query_base| **query_base == ref_base as u8)
            .count();
        let query_clean = query_repeats == qseq.len();
        Self {
            rstart,
            rend,
            qstart,
            qend,
            ref_base,
            ref_repeats,
            query_repeats,
            qseq,
            query_clean,
        }
    }
}

pub fn extract_poly_locus_info_from_record(
    record: &BamRecord,
    ref_poly_region: &Vec<(usize, usize, u8)>,
) -> Option<Vec<RefPolyLocusInfo>> {
    // if record.is_secondary() || record.is_unmapped() || record.is_supplementary() {
    //     continue;
    // }
    let record_ext = BamRecordExt::new(&record);

    let ref_start = record_ext.reference_start() as i64;
    let ref_end = record_ext.reference_end() as i64;

    let query_end = record_ext.query_alignment_end() as i64;

    let mut rpos_cursor = None;
    let mut qpos_cursor = None;

    let query_seq = record.seq().as_bytes();
    let query_str = unsafe { String::from_utf8_unchecked(query_seq.clone()) };

    let mut poly_info = vec![];

    let aligned_pairs_full: Vec<[Option<i64>; 2]> =
        record.aligned_pairs_full().into_iter().collect::<Vec<_>>();

    let mut cur_poly_region_idx =
        if let Some(poly_idx) = move_poly_idx(ref_poly_region, 0, ref_start as usize) {
            if (ref_start as usize) > ref_poly_region[poly_idx].0 {
                poly_idx + 1
            } else {
                poly_idx
            }
        } else {
            return None;
        };

    if cur_poly_region_idx >= ref_poly_region.len() {
        return None;
    }

    let mut query_poly_seq = String::new();

    for [qpos, rpos] in aligned_pairs_full.iter() {
        let qpos = *qpos;
        let rpos = *rpos;

        if qpos.is_some() {
            qpos_cursor = qpos;
        }
        if rpos.is_some() {
            rpos_cursor = rpos;
        }

        if let Some(rpos_cursor_) = rpos_cursor {
            if rpos_cursor_ < ref_start {
                continue;
            }
            if rpos_cursor_ >= ref_end {
                break;
            }
        } else {
            continue;
        }

        if let Some(qpos_cursor_) = qpos_cursor {
            if qpos_cursor_ >= query_end {
                break;
            }
        } else {
            continue;
        }

        // 收尾
        if let Some(rpos_) = rpos.map(|v| v as usize) {
            let cur_poly_region = unsafe { ref_poly_region.get_unchecked(cur_poly_region_idx) };
            if rpos_ == (cur_poly_region.1 - 1) {
                if let Some(qpos_) = qpos.map(|v| v as usize) {
                    query_poly_seq.push_str(&query_str[qpos_..qpos_ + 1]);
                }

                // poly region end
                let qend = qpos_cursor.unwrap() as usize + 1;
                let qstart = qend - query_poly_seq.len();
                poly_info.push(RefPolyLocusInfo::new(
                    cur_poly_region.0,
                    cur_poly_region.1,
                    qstart,
                    qend,
                    cur_poly_region.2 as char,
                    cur_poly_region.1 - cur_poly_region.0,
                    query_poly_seq,
                ));

                // move
                cur_poly_region_idx += 1;
                if cur_poly_region_idx >= ref_poly_region.len() {
                    break;
                }

                query_poly_seq = String::new();
            }
        }

        let cur_region = unsafe { ref_poly_region.get_unchecked(cur_poly_region_idx) };
        // region 开始 & 持续
        let rpos_cur_or_pre = rpos_cursor.unwrap() as usize;

        if let Some(rpos_) = rpos.map(|v| v as usize) {
            if rpos_ < cur_region.0 {
                continue;
            }
            // rpos_ in region
            if let Some(qpos_) = qpos.map(|v| v as usize) {
                query_poly_seq.push_str(&query_str[qpos_..qpos_ + 1]);
            }
        } else {
            if (rpos_cur_or_pre + 1) < cur_region.0 {
                continue;
            }
            if let Some(qpos_) = qpos.map(|v| v as usize) {
                query_poly_seq.push_str(&query_str[qpos_..qpos_ + 1]);
            }
        }
    }

    Some(poly_info)
}

fn move_poly_idx(
    ref_poly_region: &Vec<(usize, usize, u8)>,
    cur_idx: usize,
    cur_pos: usize,
) -> Option<usize> {
    for idx in cur_idx..ref_poly_region.len() {
        match position_relation(&ref_poly_region[cur_idx], cur_pos) {
            PosRelation::Left => return Some(idx),
            PosRelation::Middle => return Some(idx),
            PosRelation::Right => (),
        }
    }

    return None;
}

pub fn find_poly_n_regions(sequence: &[u8]) -> Vec<(usize, usize, u8)> {
    if sequence.is_empty() {
        return Vec::new();
    }

    let mut regions = Vec::new();
    let mut start = 0;
    let mut bases = sequence.iter();
    let mut pre_base = *bases.next().unwrap();

    for (mut i, &cur_base) in bases.enumerate() {
        i += 1;
        if cur_base != pre_base {
            if (i - start) > 1 {
                regions.push((start, i, pre_base));
            }
            start = i;
            pre_base = cur_base;
        }
    }

    // 保存最后一个 homopolymer 区域
    if sequence.len() - start > 1 {
        regions.push((start, sequence.len(), pre_base));
    }

    regions
}

#[derive(Debug, Clone, Copy)]
pub enum PosRelation {
    Left,
    Right,
    Middle,
}

pub fn position_relation(seb: &(usize, usize, u8), pos: usize) -> PosRelation {
    return if pos < seb.0 {
        PosRelation::Left
    } else if pos >= seb.1 {
        PosRelation::Right
    } else {
        PosRelation::Middle
    };
}

#[cfg(test)]
mod test {
    use crate::poly_n::find_poly_n_regions;

    #[test]
    fn test_find_homopolymer_regions() {
        let seq = b"AAAACCCGGTT";
        let res = find_poly_n_regions(seq);
        assert_eq!(res, vec![(0, 4, 65), (4, 7, 67), (7, 9, 71), (9, 11, 84)]);

        let seq = b"AAAACCCGGT";
        let res = find_poly_n_regions(seq);
        assert_eq!(res, vec![(0, 4, 65), (4, 7, 67), (7, 9, 71)]);

        let seq = b"ACCCGGT";
        let res = find_poly_n_regions(seq);
        assert_eq!(res, vec![(1, 4, 67), (4, 6, 71)]);

        let seq = b"AACGGT";
        let res = find_poly_n_regions(seq);
        // println!("{:?}", res);
        assert_eq!(res, vec![(0, 2, 67), (3, 5, 71)]);
    }
}
