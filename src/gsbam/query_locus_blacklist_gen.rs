use std::collections::HashSet;

use crate::itertools::sliding_window;

use super::{
    bam_record_ext::BamRecord,
    cigar_ext::{long_ins_regions_in_query, RangeIdentityCalculator},
};

pub trait TQueryLocusBlacklist: Send + Sync {
    fn get_blacklist_locus(&self, record: &BamRecord) -> HashSet<usize>;
}

/// treat the ins >= ins_thr as blacklist locus
pub struct LongInsBlacklist {
    ins_thr: usize,
}

impl LongInsBlacklist {
    pub fn new(ins_thr: usize) -> Self {
        Self { ins_thr }
    }
}

impl TQueryLocusBlacklist for LongInsBlacklist {
    fn get_blacklist_locus(&self, record: &BamRecord) -> HashSet<usize> {
        let regions = long_ins_regions_in_query(&record.cigar().take(), self.ins_thr);
        regions
            .into_iter()
            .flat_map(|(start, end)| (start..end).into_iter())
            .collect()
    }
}

pub struct LowIdentityBlacklist {
    identity_thr: f32,
    win_size: usize,
    win_ovlp: usize,
}

impl LowIdentityBlacklist {
    pub fn new(identity_thr: f32, win_size: usize, win_ovlp: usize) -> Self {
        assert!(
            win_size > win_ovlp,
            "win_size > win_ovlp, but got {} <= {}",
            win_size,
            win_ovlp
        );
        Self {
            identity_thr,
            win_size,
            win_ovlp,
        }
    }
}

impl TQueryLocusBlacklist for LowIdentityBlacklist {
    fn get_blacklist_locus(&self, record: &BamRecord) -> HashSet<usize> {
        let identity_calc = RangeIdentityCalculator::new(&record.cigar().take());
        let seq_len = record.seq_len();
        sliding_window(seq_len, self.win_size, self.win_ovlp, true)
            .into_iter()
            .flat_map(|(start, end)| {
                let (identity_start, identity_end, identity) =
                    identity_calc.compute_range_identity(start as u32, end as u32);
                // println!("{}, {}, {}", identity_start, identity_end, identity);
                
                let (identity_start, identity_end) =
                    (identity_start as usize, identity_end as usize);
                if (identity_end - identity_start) > self.win_size / 2
                    && identity < self.identity_thr
                {
                    (identity_start..identity_end).into_iter()
                } else {
                    (0..0).into_iter()
                }
            })
            .collect()
    }
}

pub fn get_query_locus_blacklist(
    record: &BamRecord,
    query_locus_blacklist_gen: Option<&Vec<Box<dyn TQueryLocusBlacklist>>>,
) -> HashSet<usize> {
    let query_locus_blacklist = query_locus_blacklist_gen
        .map(|gens| {
            gens.iter()
                .flat_map(|gen| gen.get_blacklist_locus(record).into_iter())
                .collect()
        })
        .unwrap_or(HashSet::new());
    query_locus_blacklist
}
