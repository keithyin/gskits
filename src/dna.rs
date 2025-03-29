use std::mem::MaybeUninit;

pub static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'-' as usize] = b'-';
    table[b'*' as usize] = b'*';
    table[b'N' as usize] = b'N';
    table
};

/// acgt ACGT 0 1 2 3
pub const SEQ_NT4_TABLE: [u8; 256] = [
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    let mut result: Vec<MaybeUninit<u8>> = vec![MaybeUninit::uninit(); dna.len()];
    dna.iter()
        .rev()
        .zip(result.iter_mut())
        .for_each(|(&base, slot)| {
            let complement = COMPLEMENT_TABLE[base as usize];
            slot.write(complement);
        });
    unsafe { std::mem::transmute(result) }
}

#[cfg(test)]
mod test {
    use crate::dna::reverse_complement;

    #[test]
    fn test_reverse_complement() {
        let dna_sequence = "ATCGTAGC";
        let res = reverse_complement(dna_sequence.as_bytes());
        assert_eq!(res, b"GCTACGAT");
    }

    #[test]
    fn test_bio_alpha() {
        let iupac = bio::alphabets::dna::iupac_alphabet();
        println!("{:?}", iupac);
    }
}
