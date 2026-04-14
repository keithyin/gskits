use std::{
    arch::x86_64::{
        self, __m128i, __m256i, _mm256_castps_si256, _mm256_castsi256_ps,
        _mm256_permute2x128_si256, _mm256_set1_epi16, _mm256_set_epi64x, _mm256_setr_m128i,
        _mm256_shuffle_epi32, _mm256_shuffle_epi8, _mm256_shuffle_ps, _mm256_shufflehi_epi16,
        _mm256_shufflelo_epi16, _mm256_storeu2_m128i, _mm256_unpackhi_epi16, _mm256_unpackhi_epi32,
        _mm256_unpackhi_epi64, _mm256_unpackhi_epi8, _mm256_unpacklo_epi16, _mm256_unpacklo_epi32,
        _mm256_unpacklo_epi64, _mm256_unpacklo_epi8,
    },
    mem::MaybeUninit,
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

pub fn transpose_withoutput<T>(matrix: &Vec<T>, transposed: &mut Vec<T>, dim0: usize, dim1: usize)
where
    T: Sized + Default + Clone + Copy,
{
    // 创建一个新的 Vec 用于存储转置后的矩阵

    // 遍历原矩阵并填充转置后的矩阵
    for i in 0..dim0 {
        for j in 0..dim1 {
            unsafe {
                *transposed.get_unchecked_mut(j * dim0 + i) = *matrix.get_unchecked(i * dim1 + j);
            }
        }
    }
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

pub fn mm256_shuffle_epi32<const MASK: i32>(a: __m256i, b: __m256i) -> __m256i {
    let a = unsafe { _mm256_castsi256_ps(a) };
    let b = unsafe { _mm256_castsi256_ps(b) };
    let c = unsafe { _mm256_shuffle_ps::<MASK>(a, b) };
    let c = unsafe { _mm256_castps_si256(c) };
    c
}

///
/// a: a0 a1 a2 a3 a4 a5 a6 a7
/// b: b0 b1 b2 b3 b4 b5 b6 b7
/// c: a0 b0 a2 b2 a4 b4 a6 b6
pub fn mm256_even_interleave_epi32(a: __m256i, b: __m256i) -> __m256i {
    let a = unsafe { _mm256_shuffle_epi32::<0b11011000>(a) };
    let b = unsafe { _mm256_shuffle_epi32::<0b11011000>(b) };
    unsafe { _mm256_unpacklo_epi32(a, b) }
}

///
/// a: a0 a1 a2 a3 a4 a5 a6 a7
/// b: b0 b1 b2 b3 b4 b5 b6 b7
/// c: a1 b1 a3 b3 a5 b5 a7 b7
pub fn mm256_odd_interleave_epi32(a: __m256i, b: __m256i) -> __m256i {
    let a = unsafe { _mm256_shuffle_epi32::<0b11011000>(a) };
    let b = unsafe { _mm256_shuffle_epi32::<0b11011000>(b) };
    unsafe { _mm256_unpackhi_epi32(a, b) }
}

///
/// a: a0 a1 a2 a3 a4 a5 a6 a7 ...
/// b: b0 b1 b2 b3 b4 b5 b6 b7 ...
/// c: a0 b0 a2 b2 a4 b4 a6 b6 ...
pub fn mm256_even_interleave_epi16(mut a: __m256i, mut b: __m256i) -> __m256i {
    a = unsafe {
        _mm256_shuffle_epi8(
            a,
            _mm256_set_epi64x(
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
            ),
        )
    };
    b = unsafe {
        _mm256_shuffle_epi8(
            b,
            _mm256_set_epi64x(
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
            ),
        )
    };
    unsafe { _mm256_unpacklo_epi16(a, b) }
}

///
/// a: a0 a1 a2 a3 a4 a5 a6 a7 ...
/// b: b0 b1 b2 b3 b4 b5 b6 b7 ...
/// c: a1 b1 a3 b3 a5 b5 a7 b7 ...
pub fn mm256_odd_interleave_epi16(mut a: __m256i, mut b: __m256i) -> __m256i {
    a = unsafe {
        _mm256_shuffle_epi8(
            a,
            _mm256_set_epi64x(
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
            ),
        )
    };
    b = unsafe {
        _mm256_shuffle_epi8(
            b,
            _mm256_set_epi64x(
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
                0x0f0e0b0a07060302,
                0x0d0c090805040100,
            ),
        )
    };
    unsafe { _mm256_unpackhi_epi16(a, b) }
}

/// a: a0 a1 a2 a3 a4 a5 a6 a7 ...
/// b: b0 b1 b2 b3 b4 b5 b6 b7 ...
/// c: a0 b0 a2 b2 a4 b4 a6 b6 ...
pub fn mm256_even_interleave_epi8(mut a: __m256i, mut b: __m256i) -> __m256i {
    a = unsafe {
        _mm256_shuffle_epi8(
            a,
            _mm256_set_epi64x(
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
            ),
        )
    };
    b = unsafe {
        _mm256_shuffle_epi8(
            b,
            _mm256_set_epi64x(
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
            ),
        )
    };
    unsafe { _mm256_unpacklo_epi8(a, b) }
}

///
/// a: a0 a1 a2 a3 a4 a5 a6 a7 ...
/// b: b0 b1 b2 b3 b4 b5 b6 b7 ...
/// c: a1 b1 a3 b3 a5 b5 a7 b7 ...
pub fn mm256_odd_interleave_epi8(mut a: __m256i, mut b: __m256i) -> __m256i {
    a = unsafe {
        _mm256_shuffle_epi8(
            a,
            _mm256_set_epi64x(
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
            ),
        )
    };
    b = unsafe {
        _mm256_shuffle_epi8(
            b,
            _mm256_set_epi64x(
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
                0x0f0d0b0907050301,
                0x0e0c0a0806040200,
            ),
        )
    };
    unsafe { _mm256_unpackhi_epi8(a, b) }
}

fn avx2_transopose_u8_16x16_64(rows: [__m256i; 8]) -> [__m256i; 8] {
    let mut result_rows = [MaybeUninit::uninit(); 8];
    result_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        if idx < 4 {
            v.write(unsafe { _mm256_unpacklo_epi64(rows[idx], rows[idx + 4]) });
        } else {
            v.write(unsafe { _mm256_unpackhi_epi64(rows[idx - 4], rows[idx]) });
        }
    });

    unsafe { std::mem::transmute(result_rows) }
}

fn avx2_transopose_u8_16x16_32(rows: [__m256i; 8]) -> [__m256i; 8] {
    let mut result_rows = [MaybeUninit::uninit(); 8];
    result_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        if idx % 4 < 2 {
            v.write(mm256_even_interleave_epi32(rows[idx], rows[idx + 2]));
        } else {
            v.write(mm256_odd_interleave_epi32(rows[idx - 2], rows[idx]));
        }
    });
    unsafe { std::mem::transmute(result_rows) }
}

fn avx2_transopose_u8_16x16_16(rows: [__m256i; 8]) -> [__m256i; 8] {
    let mut result_rows = [MaybeUninit::uninit(); 8];
    result_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        if idx % 2 == 0 {
            v.write(mm256_even_interleave_epi16(rows[idx], rows[idx + 1]));
        } else {
            v.write(mm256_odd_interleave_epi16(rows[idx - 1], rows[idx]));
        }
    });
    unsafe { std::mem::transmute(result_rows) }
}

fn avx2_transopose_u8_16x16_8(rows: [__m256i; 8]) -> [__m256i; 8] {
    let mut result_rows = [MaybeUninit::uninit(); 8];
    result_rows.iter_mut().enumerate().for_each(|(idx, v)| {
        let hi_2_lo = unsafe { _mm256_permute2x128_si256::<0x01>(rows[idx], rows[idx]) };
        let tmp1 = mm256_even_interleave_epi8(rows[idx], hi_2_lo);
        let tmp2 = mm256_odd_interleave_epi8(rows[idx], hi_2_lo);
        let res = unsafe { _mm256_permute2x128_si256::<0b00100000>(tmp1, tmp2) };
        v.write(res);
    });

    unsafe { std::mem::transmute(result_rows) }
}

fn avx_transpose_u8_store(rows: [__m256i; 8], output: &mut [&mut [u8]]) {
    rows.into_iter().enumerate().for_each(|(idx, v)| unsafe {
        _mm256_storeu2_m128i(
            output[idx * 2 + 1].as_mut_ptr() as *mut __m128i,
            output[idx * 2].as_mut_ptr() as *mut __m128i,
            rows[idx],
        );
    });
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

    let rows: [__m256i; 8] = unsafe { std::mem::transmute(origin_rows) };

    let rows = avx2_transopose_u8_16x16_64(rows);
    // avx_transpose_u8_store(rows, output);

    // println!("out:\n{:?}", output);
    let rows = avx2_transopose_u8_16x16_32(rows);
    // avx_transpose_u8_store(rows, output);

    // println!("out:\n{:?}", output);
    let rows = avx2_transopose_u8_16x16_16(rows);
    // avx_transpose_u8_store(rows, output);

    // println!("out:\n{:?}", output);
    let rows = avx2_transopose_u8_16x16_8(rows);
    avx_transpose_u8_store(rows, output);

    // println!("out:\n{:?}", output);
}

#[cfg(test)]
mod test {
    use std::arch::x86_64::{__m256i, _mm256_loadu_si256, _mm256_storeu_si256};

    use crate::matrix::transpose::{
        avx2_transopose_u8_16x16, mm256_even_interleave_epi16, mm256_even_interleave_epi32,
        mm256_odd_interleave_epi16, mm256_odd_interleave_epi32,
    };

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

    #[test]
    fn test_even_interleave() {
        unsafe {
            // 输入数据：a = [0,1,2,3,4,5,6,7], b = [10,11,12,13,14,15,16,17]
            let a_vals = [0i32, 1, 2, 3, 4, 5, 6, 7];
            let b_vals = [10i32, 11, 12, 13, 14, 15, 16, 17];
            let a = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
            let b = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);

            let result = mm256_even_interleave_epi32(a, b);

            let mut out = [0i32; 8];
            _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, result);

            let expected = [0, 10, 2, 12, 4, 14, 6, 16];
            assert_eq!(out, expected);
        }
    }

    #[test]
    fn test_odd_interleave() {
        unsafe {
            // 输入数据：a = [0,1,2,3,4,5,6,7], b = [10,11,12,13,14,15,16,17]
            let a_vals = [0i32, 1, 2, 3, 4, 5, 6, 7];
            let b_vals = [10i32, 11, 12, 13, 14, 15, 16, 17];
            let a = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
            let b = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);

            let result = mm256_odd_interleave_epi32(a, b);

            let mut out = [0i32; 8];
            _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, result);

            let expected = [1, 11, 3, 13, 5, 15, 7, 17];
            assert_eq!(out, expected);
        }
    }

    #[test]
    fn test_even_interleave_epi16() {
        unsafe {
            // 输入 16 位整数：a = 0,1,2,...,15；b = 100,101,102,...,115
            let mut a_vals = [0i16; 16];
            let mut b_vals = [0i16; 16];
            for i in 0..16 {
                a_vals[i] = i as i16;
                b_vals[i] = 100 + i as i16;
            }
            let a = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
            let b = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);

            let result = mm256_even_interleave_epi16(a, b);

            let mut out = [0i16; 16];
            _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, result);

            // 预期：偶数索引取自 a，然后来自 b 的相同位置（偶数索引）
            let mut expected = [0i16; 16];
            for i in 0..8 {
                expected[2 * i] = a_vals[2 * i];
                expected[2 * i + 1] = b_vals[2 * i];
            }
            assert_eq!(out, expected);
        }
    }

    #[test]
    fn test_odd_interleave_epi16() {
        unsafe {
            // 输入 16 位整数：a = 0,1,2,...,15；b = 100,101,102,...,115
            let mut a_vals = [0i16; 16];
            let mut b_vals = [0i16; 16];
            for i in 0..16 {
                a_vals[i] = i as i16;
                b_vals[i] = 100 + i as i16;
            }
            let a = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
            let b = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);

            let result = mm256_odd_interleave_epi16(a, b);

            let mut out = [0i16; 16];
            _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, result);

            // 预期：奇数索引取自 a，然后来自 b 的相同奇数索引
            let mut expected = [0i16; 16];
            for i in 0..8 {
                expected[2 * i] = a_vals[2 * i + 1];
                expected[2 * i + 1] = b_vals[2 * i + 1];
            }
            assert_eq!(out, expected);
        }
    }

    #[test]
    fn test_mm256_even_interleave_epi8() {
        // 固定测试：使用递增序列，便于人工验证
        let mut a_vals = [0u8; 32];
        let mut b_vals = [0u8; 32];
        for i in 0..32 {
            a_vals[i] = i as u8;
            b_vals[i] = (i + 32) as u8;
        }
        let a = unsafe { _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i) };
        let b = unsafe { _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i) };
        let c = super::mm256_even_interleave_epi8(a, b);

        let mut result = [0u8; 32];
        unsafe { _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, c) };

        // 预期结果：a0,b0,a2,b2,...,a30,b30
        let mut expected = [0u8; 32];
        for i in 0..16 {
            expected[2 * i] = a_vals[2 * i];
            expected[2 * i + 1] = b_vals[2 * i];
        }
        assert_eq!(result, expected, "固定模式测试失败");

        // 随机测试：多组随机数据
        // 使用简单的 xorshift 生成随机数，避免外部依赖
        let mut rng_state = 123456789u64;
        let mut rand_byte = || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            (rng_state & 0xFF) as u8
        };

        for _ in 0..100 {
            // 生成随机 a 和 b
            for i in 0..32 {
                a_vals[i] = rand_byte();
                b_vals[i] = rand_byte();
            }
            let a = unsafe { _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i) };
            let b = unsafe { _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i) };
            let c = super::mm256_even_interleave_epi8(a, b);
            unsafe { _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, c) };

            // 计算预期
            for i in 0..16 {
                expected[2 * i] = a_vals[2 * i];
                expected[2 * i + 1] = b_vals[2 * i];
            }
            assert_eq!(result, expected, "随机测试失败");
        }
    }

    fn to_m256i(arr: &[u8; 32]) -> __m256i {
        unsafe { std::arch::x86_64::_mm256_loadu_si256(arr.as_ptr() as *const __m256i) }
    }

    /// 安全地将 __m256i 转换回 32 字节数组。
    fn from_m256i(vec: __m256i) -> [u8; 32] {
        let mut arr = [0u8; 32];
        unsafe { std::arch::x86_64::_mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, vec) };
        arr
    }

    #[test]
    fn test_mm256_odd_interleave_epi8() {
        // 固定测试：使用递增序列
        let a_vals: [u8; 32] = core::array::from_fn(|i| i as u8);
        let b_vals: [u8; 32] = core::array::from_fn(|i| (i + 32) as u8);

        let a = to_m256i(&a_vals);
        let b = to_m256i(&b_vals);
        let c = super::mm256_odd_interleave_epi8(a, b);
        let result = from_m256i(c);

        // 预期结果：a1, b1, a3, b3, ..., a31, b31
        let mut expected = [0u8; 32];
        for i in 0..16 {
            expected[2 * i] = a_vals[2 * i + 1];
            expected[2 * i + 1] = b_vals[2 * i + 1];
        }
        assert_eq!(result, expected, "固定模式测试失败");

        // 随机测试：多组随机数据
        let mut rng_state = 123456789u64;
        let mut rand_byte = || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            (rng_state & 0xFF) as u8
        };

        for _ in 0..100 {
            let mut a_rand = [0u8; 32];
            let mut b_rand = [0u8; 32];
            for i in 0..32 {
                a_rand[i] = rand_byte();
                b_rand[i] = rand_byte();
            }
            let a = to_m256i(&a_rand);
            let b = to_m256i(&b_rand);
            let c = super::mm256_odd_interleave_epi8(a, b);
            let result = from_m256i(c);

            let mut expected = [0u8; 32];
            for i in 0..16 {
                expected[2 * i] = a_rand[2 * i + 1];
                expected[2 * i + 1] = b_rand[2 * i + 1];
            }
            assert_eq!(result, expected, "随机测试失败");
        }
    }
}
