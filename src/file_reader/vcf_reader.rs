use std::{collections::HashMap, io::{self, BufRead}};

#[derive(Debug)]
pub struct VcfRowData {
    pub chrom: String,
    pub pos: usize, // zero based. 1-based in vcf file
    pub ref_bases: String,
    pub alt_bases: String,
    pub phred_q: u32,
}

impl VcfRowData {

    /// chrom pos id ref alt qual filter info format ..
    pub fn from_str(inp: &str) -> Self {
        let items = inp.split("\t").collect::<Vec<_>>();

        let chrom = items[0].trim().to_string();
        let pos = items[1].parse::<usize>().unwrap() - 1;
        let ref_bases = items[3].trim().to_string();
        let alt_bases = items[4].trim().to_string();
        let phred_q = items[5].parse::<u32>().unwrap();

        VcfRowData {
            chrom,
            pos, 
            ref_bases, 
            alt_bases, 
            phred_q
        }
    }

    #[allow(unused)]
    pub fn new(chrom: String, pos: usize, ref_bases: String, alt_bases: String, phred_q: u32) -> Self {
        Self { 
            chrom, 
            pos, 
            ref_bases, 
            alt_bases, 
            phred_q
        }
    }

}

pub struct VcfReaderIter<'a> {
    reader: &'a mut dyn io::BufRead
}

impl <'a> VcfReaderIter<'a>{
    
    pub fn new(reader: &'a mut dyn io::BufRead) -> Self {
        VcfReaderIter { reader: reader }
    }

}

impl<'a> Iterator for VcfReaderIter<'a> {
    type Item = VcfRowData;

    fn next(&mut self) -> Option<Self::Item> {
        
        let mut line = String::new();
        loop {
            line.clear();
            if let Ok(n) = self.reader.read_line(&mut line) {
                if n == 0 {
                    return None;
                }
                if line.starts_with("#") {
                    continue;
                }

                return Some(VcfRowData::from_str(&line));
            } else {
                return None;
            }

        }

    }

}


pub struct VcfInfo {
    info: HashMap<String, Vec<usize>>,
}

impl VcfInfo {
    pub fn new(file_reader: &mut dyn BufRead) -> Self {
        let mut info = HashMap::new();
        let file_iter = VcfReaderIter::new(file_reader);
        file_iter.for_each(|vcf| {
            // .into_iter().collect::<Vec<usize>>()
            let bad_points = vcf.pos..(vcf.pos + vcf.ref_bases.len()) ;
            info.entry(vcf.chrom).or_insert(vec![]).extend(bad_points);

        });

        info.iter_mut().for_each(|(_k, v)| {
            v.sort_unstable();
        });

        Self { 
            info
        }
    }

    /// true: means the range contains a variant loci
    pub fn range_hit(&self, chromosome: &str, range: &(usize, usize)) -> bool {
        
        if !self.info.contains_key(chromosome) {
            return false;
        }

        let sorted_vec = self.info.get(chromosome).unwrap();
        return match sorted_vec.binary_search(&range.0) {
            Ok(_) => true,
            Err(n) => {
                if n == sorted_vec.len() {
                    false
                } else {
                    if range.1 <= sorted_vec[n] {
                        false
                    } else {
                        true
                    }
                }
            }
        };
    }

}

