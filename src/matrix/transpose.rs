use std::{
    arch::x86_64::{self, __m128i, __m256i, _mm256_storeu2_m128i},
    mem::{ManuallyDrop, MaybeUninit},
};

pub fn transpose<T>(matrix: &Vec<T>, dim0: usize, dim1: usize) -> Vec<T>
where
    T: Sized + Default + Clone + Copy,
{
    // 创建一个新的 Vec 用于存储转置后的矩阵
    let mut transposed = vec![T::default(); dim0 * dim1];

    // 遍历原矩阵并填充转置后的矩阵
    for i in 0..dim0 {
        for j in 0..dim1 {
            unsafe {
                *transposed.get_unchecked_mut(j * dim0 + i) = *matrix.get_unchecked(i * dim1 + j);
            }
        }
    }

    transposed
}

pub fn transpose_blocked<T>(
    matrix: &Vec<T>,
    dim0: usize,
    dim1: usize,
    block_size: Option<usize>,
) -> Vec<T>
where
    T: Sized + Default + Clone + Copy,
{
    assert_eq!(
        matrix.len(),
        dim0 * dim1,
        "matrix.len={}, dim0*dim1={}",
        matrix.len(),
        dim0 * dim1
    );

    let block_size = block_size.unwrap_or(500);
    // 创建一个新的 Vec 用于存储转置后的矩阵
    let mut transposed = vec![T::default(); dim0 * dim1];

    // 分块转置
    for i_block in (0..dim0).step_by(block_size) {
        for j_block in (0..dim1).step_by(block_size) {
            // 转置每个小块
            let i_end = (i_block + block_size).min(dim0);
            let j_end = (j_block + block_size).min(dim1);

            for i in i_block..i_end {
                for j in j_block..j_end {
                    unsafe {
                        *transposed.get_unchecked_mut(j * dim0 + i) =
                            *matrix.get_unchecked(i * dim1 + j);
                    }
                }
            }
        }
    }

    transposed
}

pub fn avx2_transopose_u8_16x16(input: &[&[u8]], output: &mut [&mut [u8]]) {
    let mut origin_rows = [MaybeUninit::uninit(); 8];
    origin_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        v.write(unsafe {
            x86_64::_mm256_loadu2_m128i(
                input[idx * 2 + 1].as_ptr() as *const __m128i,
                input[idx * 2].as_ptr() as *const __m128i,
            )
        });
    });

    let origin_rows: [__m256i; 8] = unsafe { std::mem::transmute(origin_rows) };

    let mut tmp1_rows = [MaybeUninit::uninit(); 4];
    tmp1_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        v.write(unsafe { x86_64::_mm256_unpacklo_epi64(origin_rows[idx], origin_rows[idx + 4]) });
    });
    let tmp1_rows: [__m256i; 4] = unsafe { std::mem::transmute(tmp1_rows) };

    let mut tmp2_rows = [MaybeUninit::uninit(); 2];
    tmp2_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        v.write(unsafe { x86_64::_mm256_unpacklo_epi32(tmp1_rows[idx], origin_rows[idx + 2]) });
    });
    let tmp2_rows: [__m256i; 2] = unsafe { std::mem::transmute(tmp2_rows) };

    let tmp3_row = unsafe { x86_64::_mm256_unpacklo_epi16(tmp2_rows[0], tmp2_rows[1]) };
    let hi = unsafe {
        x86_64::_mm256_permute2x128_si256::<0x01>(tmp3_row, tmp3_row)
    };
    let ok_res = unsafe {
        x86_64::_mm256_unpacklo_epi8(tmp3_row, hi)
    };

    

    // let mut tmp2_rows =

    unsafe {
        _mm256_storeu2_m128i(
            output[1].as_mut_ptr() as *mut __m128i,
            output[0].as_mut_ptr() as *mut __m128i,
            tmp3_row,
        );

        _mm256_storeu2_m128i(
            output[3].as_mut_ptr() as *mut __m128i,
            output[2].as_mut_ptr() as *mut __m128i,
            hi,
        );

        _mm256_storeu2_m128i(
            output[5].as_mut_ptr() as *mut __m128i,
            output[4].as_mut_ptr() as *mut __m128i,
            ok_res,
        );
    }

    println!("out:{:?}", output);
}

#[cfg(test)]
mod test {
    use crate::matrix::transpose::avx2_transopose_u8_16x16;

    use super::{transpose, transpose_blocked};

    #[test]
    fn test_transpose() {
        let values = (0..8).into_iter().collect::<Vec<_>>();
        let res = transpose(&values, 2, 4);
        assert_eq!(res, vec![0, 4, 1, 5, 2, 6, 3, 7]);
    }

    #[test]
    fn test_transpose_blocked() {
        let values = (0..8).into_iter().collect::<Vec<_>>();
        let res1 = transpose(&values, 2, 4);
        let res2 = transpose_blocked(&values, 2, 4, Some(3));
        assert_eq!(res1, res2);
    }

    #[test]
    fn test_transpose_avx2() {
        let value = (0..(16 * 16))
            .into_iter()
            .map(|v| v as u8)
            .collect::<Vec<u8>>();
        let value2 = (0..16)
            .into_iter()
            .map(|idx| &value[idx * 16..(idx + 1) * 16])
            .collect::<Vec<&[u8]>>();
        let mut output = (0..16)
            .into_iter()
            .map(|v| vec![0_u8; 16])
            .collect::<Vec<_>>();
        let mut output2 = output
            .iter_mut()
            .map(|v| v.as_mut_slice())
            .collect::<Vec<_>>();

        avx2_transopose_u8_16x16(&value2, &mut output2);
    }
}
