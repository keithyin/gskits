use std::{collections::HashMap, io::BufRead};


#[derive(Debug, Clone)]
pub struct BedRowData {
    pub chromosome: String,
    pub begin: usize, // 0-coordinate
    pub end: usize
}

impl BedRowData {
    pub fn from_str(inp: &str) -> Self {
        let items = inp.trim().split("\t").collect::<Vec<_>>();
        let chromosome = items[0].to_string();
        let begin = items[1].parse::<usize>().unwrap();
        let end = items[2].parse::<usize>().unwrap();
        Self { 
            chromosome, 
            begin,
            end 
        }
    }
}

pub struct BedFileReaderIter<'a> {
    reader: &'a mut dyn BufRead
}

impl<'a> BedFileReaderIter<'a>{
    pub fn new(reader: &'a mut dyn BufRead) -> Self {
        Self { 
            reader: reader 
        }
    }
}

impl<'a> Iterator for BedFileReaderIter<'a>  {
    type Item = BedRowData;
    fn next(&mut self) -> Option<Self::Item> {
        
        let mut line = String::new();
        if let Ok(n) = self.reader.read_line(&mut line) {
            if n == 0 {
                return None;
            }

            return Some(BedRowData::from_str(&line));

        } else {
            return None;
        }

    }
}

///
/// bed file
///     chrom   start   end   [0-based]
pub struct BedInfo {
    regions: HashMap<String, Vec<(usize, usize)>>
}

impl BedInfo {
    pub fn new(file_reader: &mut dyn BufRead) -> Self {
        let mut info = HashMap::new();
        let reader_iter = BedFileReaderIter::new(file_reader);
        reader_iter.for_each(|bed_row| {
            info.entry(bed_row.chromosome)
                .or_insert(vec![])
                .push((bed_row.begin, bed_row.end));
        });

        info.iter_mut().for_each(|(_k, v)| {
            v.sort_unstable_by_key(|v| v.0);
        });

        Self { regions: info}
    }

    #[allow(unused)]
    pub fn from_info(info: HashMap<String, Vec<(usize, usize)> >) -> Self {
        Self { regions: info }
    }

    pub fn within_the_range(&self, chromosome: &str, begin_end: &(usize, usize)) -> bool{
        
        if !self.regions.contains_key(chromosome) {
            return false;
        }
        let sorted_range = self.regions.get(chromosome).unwrap();
        return match sorted_range.binary_search_by_key(&begin_end.0, |v| v.0) {
            Ok(n) => sorted_range[n].1 >= begin_end.1,
            Err(n) => {
                if n == sorted_range.len() || n == 0{
                    false
                } else {
                    begin_end.1 <= sorted_range[n-1].1
                }
            }
        };
    }

    pub fn point_within_region(&self, chrom: &str, position: usize) -> bool {
        self.regions.get(chrom)
        .and_then(|start_ends| {
            let in_region = match start_ends.binary_search_by_key(&position, |&(s, _)| s) {
                Ok(_) => true,
                Err(n) => {
                    if n == 0 {
                        false
                    } else {
                        let pre_s_e = &start_ends[n - 1];
                        if position >= pre_s_e.1 {
                            false
                        } else {
                            true
                        }
                    }
                }
            };
            Some(in_region)

        })
        .unwrap_or(false)

    }

    pub fn get_regions(&self, chrom: &str) -> Option<&Vec<(usize, usize)>> {
        self.regions.get(chrom)

    }

}