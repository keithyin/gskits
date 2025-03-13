use crate::ds::ReadInfo;

pub mod fasta_reader;
pub mod fastq_reader;
pub mod fastx2bam;


pub fn fastx_header_line_to_header(header_line: &str) -> String {
    let header = if header_line[1..].contains(" ") {
        header_line[1..].split_once(" ").unwrap().0.to_string()
    } else {
        header_line[1..].to_string()
    };
    return header;
}


pub fn read_fastx<T: Iterator<Item = ReadInfo>>(reader: T) -> Vec<ReadInfo> {
    reader.into_iter().collect::<_>()
}