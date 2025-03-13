
pub mod gsbam;
pub mod fastx_reader;
pub mod pbar;
pub mod samtools;
pub mod utils;
pub mod ds;
pub mod dna;
pub mod file_reader;
pub mod phreq;
pub mod matrix;
pub mod itertools;
pub mod poly_n;
pub mod cleanup;
pub mod collections_ext;

#[cfg(test)]
mod test {

    #[test]
    fn test_alphabet() {

        use bio::alphabets;
        println!("{:?}", alphabets::dna::alphabet());

    }

}