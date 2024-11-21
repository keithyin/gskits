use rust_htslib::bam::{self, Record};

use crate::utils::generate_tmp_filename;
use crate::fastx_reader::{fasta_reader, fastq_reader};

pub fn fasta2bam(fa_filename: &str, delim: &str, channel_idx: usize) -> String {

    let fasta_reader = fasta_reader::FastaFileReader::new(fa_filename.to_string());

    let mut header = bam::Header::new();
    let mut hd = bam::header::HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");
    header.push_record(&hd);

    let mut hd = bam::header::HeaderRecord::new(b"SQ");
    hd.push_tag(b"SN", "chr1");
    hd.push_tag(b"LN", "1234");
    header.push_record(&hd);


    let o_filepath = format!("{}.bam", fa_filename.rsplit_once(".").unwrap().0);
    let o_filepath = generate_tmp_filename(&o_filepath);
    let mut bam_writer = bam::Writer::from_path(&o_filepath, &header, bam::Format::Bam).unwrap();
    bam_writer.set_threads(4).unwrap();

    for read_info in fasta_reader {
        let qname = read_info.name;
        let seq = read_info.sequence;
        let items = qname.split(delim).skip(channel_idx).collect::<Vec<_>>();
        let ch = items[0].parse::<usize>().unwrap();

        let mut bam_record = Record::new();
        bam_record.set(qname.as_bytes(), None, seq.as_bytes(), &vec![255; seq.len()]);
        bam_record.push_aux(b"ch", bam::record::Aux::U32(ch as u32)).unwrap();
        bam_record.set_tid(0);
        bam_writer.write(&bam_record).unwrap();
    }

    o_filepath
}


pub fn fastq2bam(fq_filename: &str, delim: &str, channel_idx: usize) -> String {

    let fq_reader = fastq_reader::FastqReader::new(fq_filename.to_string());

    let mut header = bam::Header::new();
    let mut hd = bam::header::HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");
    header.push_record(&hd);

    let mut hd = bam::header::HeaderRecord::new(b"SQ");
    hd.push_tag(b"SN", "chr1");
    hd.push_tag(b"LN", "1234");
    header.push_record(&hd);


    let o_filepath = format!("{}.bam", fq_filename.rsplit_once(".").unwrap().0);
    let o_filepath = generate_tmp_filename(&o_filepath);
    let mut bam_writer = bam::Writer::from_path(&o_filepath, &header, bam::Format::Bam).unwrap();
    bam_writer.set_threads(4).unwrap();

    for read_info in fq_reader {
        let qname = read_info.name;
        let seq = read_info.sequence;
        let qual = read_info.qual.unwrap();

        let items = qname.split(delim).skip(channel_idx).collect::<Vec<_>>();
        let ch = items[0].parse::<usize>().unwrap();

        let mut bam_record = Record::new();

        let qual = qual.into_bytes().into_iter().map(|q| q-33).collect::<Vec<_>>();
        bam_record.set_tid(0);

        bam_record.set(qname.as_bytes(), None, seq.as_bytes(), &qual);
        bam_record.push_aux(b"ch", bam::record::Aux::U32(ch as u32)).unwrap();
        bam_writer.write(&bam_record).unwrap();
    }

    o_filepath
}