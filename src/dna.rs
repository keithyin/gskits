pub static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table
};

pub fn reverse_complement(dna: &str) -> String {
    let mut result = dna.as_bytes().to_vec();
    let len = result.len();

    for (i, &base) in dna.as_bytes().iter().enumerate() {
        let complement = COMPLEMENT_TABLE[base as usize];
        result[len - i - 1] = complement;
    }

    unsafe { String::from_utf8_unchecked(result) }
}

#[cfg(test)]
mod test {
    use crate::dna::reverse_complement;

    #[test]
    fn test_reverse_complement() {
        let dna_sequence = "ATCGTAGC";
        let res = reverse_complement(dna_sequence);
        assert_eq!(res, "GCTACGAT");
    }
}