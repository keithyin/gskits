use std::ops::{Deref, DerefMut};

#[derive(Debug, Clone)]
pub struct Region {
    start: usize,
    end: usize,
}

impl Region {
    pub fn new(s: usize, e: usize) -> Self {
        Self { start: s, end: e }
    }
    pub fn length(&self) -> usize {
        self.end - self.start
    }
}

pub struct Regions(Vec<Region>);

impl Deref for Regions {
    type Target = Vec<Region>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Regions {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<&Vec<(usize, usize)>> for Regions {
    fn from(value: &Vec<(usize, usize)>) -> Self {
        let regions = value
            .iter()
            .map(|v| Region {
                start: v.0,
                end: v.1,
            })
            .collect::<Vec<_>>();
        Self(regions)
    }
}

impl Regions {
    pub fn new(regions: Vec<Region>) -> Self {
        Self(regions)
    }

    pub fn total_length(&self) -> usize {
        self.iter().map(|v| v.length()).sum()
    }

    pub fn merge_regions(&self) -> Self {
        let mut regions = self.0.clone();
        regions.sort_by_key(|v| v.start);
        let mut merged_regions: Vec<Region> = Vec::new();

        for region in regions {
            if let Some(last) = merged_regions.last_mut() {
                if region.start <= last.end {
                    last.end = last.end.max(region.end);
                } else {
                    merged_regions.push(region);
                }
            } else {
                merged_regions.push(region);
            }
        }
        Self(merged_regions)
    }

    pub fn ovlp_length(&self) -> usize {
        let mut total_overlap = 0;
        let mut events = Vec::new();

        for region in &self.0 {
            events.push((region.start, 1));
            events.push((region.end, -1));
        }

        events.sort();

        let mut current_depth = 0;
        let mut prev_pos = events[0].0;

        for (pos, change) in events {
            if current_depth > 1 {
                total_overlap += pos - prev_pos;
            }
            current_depth += change;
            prev_pos = pos;
        }

        total_overlap
    }

    pub fn ovlp_ratio(&self) -> f32 {
        let merged_regions = self.merge_regions();
        let len_after_merge = merged_regions.total_length();
        let ovlp_len = self.ovlp_length();
        if len_after_merge == 0 {
            0.0
        } else {
            ovlp_len as f32 / len_after_merge as f32
        }
    }

    /// the region must be sorted for use this function
    pub fn gaps(&self, init_pos: Option<usize>, end_pos: Option<usize>) -> Vec<i64> {
        if self.is_empty() {
            return vec![];
        }
        let mut gaps = vec![];

        let mut pre_pos = init_pos;
        self.iter().for_each(|region| {
            if let Some(pre_pos) = pre_pos {
                let gap = region.start as i64 - pre_pos as i64;
                gaps.push(gap);
            }
            pre_pos = Some(region.end);
        });

        if let Some(end_pos) = end_pos {
            let gap = end_pos as i64 - pre_pos.unwrap() as i64;
            gaps.push(gap);
        }

        gaps
    }
}

#[cfg(test)]
mod test {
    use crate::ds::region::Regions;


    #[test]
    fn test_calcute_overlap_metrics() {
        let regions = vec![(1_usize, 5_usize), (1, 5)];
        let regions: Regions = (&regions).into();
        let res = regions.ovlp_ratio();
        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (5, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.ovlp_ratio();
        println!("{:?}", res);
    }

    #[test]
    fn test_calcute_gaps() {
        let regions = vec![(1_usize, 5_usize), (1, 5)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(Some(0), Some(100));
        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (5, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(Some(0), Some(100));

        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (5, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(None, None);

        println!("{:?}", res);

        let regions = vec![(1_usize, 5_usize), (7, 10)];
        let regions: Regions = (&regions).into();
        let res = regions.gaps(None, None);

        println!("{:?}", res);
    }
}
