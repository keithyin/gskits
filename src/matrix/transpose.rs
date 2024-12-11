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
    assert_eq!(matrix.len(), dim0 * dim1, "matrix.len={}, dim0*dim1={}", matrix.len(), dim0 * dim1);

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

#[cfg(test)]
mod test {
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
}
