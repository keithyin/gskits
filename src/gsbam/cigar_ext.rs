use rust_htslib::bam::record::{Cigar, CigarString};

#[derive(Debug, Clone)]
pub struct RangeIdentityCalculator {
    query_start_pos_and_op: Vec<(u32, Cigar)>,
}

impl RangeIdentityCalculator {
    pub fn new(cigar_str: &CigarString) -> Self {
        let mut q_start_cusor = 0;

        let query_start_pos_and_op = cigar_str
            .iter()
            .map(|cigar| {
                let ret = (q_start_cusor, *cigar);
                match *cigar {
                    Cigar::Match(n)
                    | Cigar::Diff(n)
                    | Cigar::SoftClip(n)
                    | Cigar::Ins(n)
                    | Cigar::Equal(n) => q_start_cusor += n,
                    _ => {}
                };
                ret
            })
            .collect::<Vec<_>>();

        Self {
            query_start_pos_and_op,
        }
    }

    pub fn range_identity(&self, start: u32, end: u32) -> f32 {
        let start_idx = match self
            .query_start_pos_and_op
            .binary_search_by_key(&start, |&(a, _)| a)
        {
            Ok(mut idx) => {
                while (idx + 1) < self.query_start_pos_and_op.len()
                    && self.query_start_pos_and_op[idx + 1].0 == start
                {
                    idx += 1;
                }
                idx
            }
            Err(idx) => idx - 1,
        };

        let end_idx = match self
            .query_start_pos_and_op
            .binary_search_by_key(&end, |&(a, _)| a)
        {
            Ok(mut idx) => {
                unsafe {
                    while idx > 0 && self.query_start_pos_and_op.get_unchecked(idx - 1).0 == end {
                        idx -= 1;
                    }
                }

                idx
            }
            Err(idx) => idx - 1,
        };

        println!("{}-{}", start_idx, end_idx);

        let res_cigars = if start_idx == end_idx {
            let res_cigar = match self.query_start_pos_and_op[start_idx].1 {
                Cigar::Diff(_) => Cigar::Diff(end - start),
                Cigar::Ins(_) => Cigar::Diff(end - start),
                Cigar::Equal(_) => Cigar::Equal(end - start),
                a => panic!("not a valid cigar in here {}", a),
            };
            vec![res_cigar]
        } else {
            let init_pos = self.query_start_pos_and_op[start_idx].0;
            let first_cigar = match self.query_start_pos_and_op[start_idx].1 {
                Cigar::Diff(n) => Cigar::Diff(n - (start - init_pos)),
                Cigar::Ins(n) => Cigar::Diff(n - (start - init_pos)),
                Cigar::Equal(n) => Cigar::Equal(n - (start - init_pos)),
                a => panic!("not a valid cigar in here {}", a),
            };

            let init_pos = self.query_start_pos_and_op[end_idx].0;
            let last_cigar = match self.query_start_pos_and_op[end_idx].1 {
                Cigar::Diff(_) => Cigar::Diff(end - init_pos),
                Cigar::Ins(_) => Cigar::Diff(end - init_pos),
                Cigar::Equal(_) => Cigar::Equal(end - init_pos),
                a => panic!("not a valid cigar in here {}", a),
            };

            let mut res_cigars = vec![first_cigar];
            res_cigars.extend(
                self.query_start_pos_and_op[start_idx + 1..end_idx]
                    .iter()
                    .map(|&(_, cigar)| cigar),
            );
            res_cigars.push(last_cigar);
            res_cigars
        };

        let span_len = res_cigars
            .iter()
            .map(|&cigar| match cigar {
                Cigar::Diff(n) | Cigar::Equal(n) | Cigar::Ins(n) => n,
                _ => 0,
            })
            .sum::<u32>();
        let eq_len = res_cigars
            .iter()
            .map(|&cigar| match cigar {
                Cigar::Equal(n) => n,
                _ => 0,
            })
            .sum::<u32>();

        let span_len = if span_len == 0 {
            1.0_f32
        } else {
            span_len as f32
        };
        eq_len as f32 / span_len
    }
}

pub fn parse_cigar_string(cigar: &str) -> Result<CigarString, String> {
    let mut cigar_ops = Vec::with_capacity(cigar.len() / 2); // 预分配容量
    let mut length = 0u32; // 累积长度

    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            // 使用 ASCII 操作来解析数字，避免字符串操作
            length = length * 10 + (c - b'0') as u32;
        } else {
            // non digit to CigarOp
            let op = match c {
                b'M' => Cigar::Match(length),
                b'=' => Cigar::Equal(length),
                b'X' => Cigar::Diff(length),
                b'I' => Cigar::Ins(length),
                b'D' => Cigar::Del(length),
                b'N' => Cigar::RefSkip(length),
                b'S' => Cigar::SoftClip(length),
                b'H' => Cigar::HardClip(length),
                b'P' => Cigar::Pad(length),
                _ => return Err(format!("Invalid CIGAR operator: {}", c as char)),
            };
            cigar_ops.push(op);
            length = 0; // reset length
        }
    }

    if length != 0 {
        return Err(format!("Trailing digits in CIGAR string: {}", cigar));
    }

    Ok(CigarString(cigar_ops))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_range_identity() {
        // 测试数据
        let cigar_str = vec![
            Cigar::Equal(5),
            Cigar::Diff(3),
            Cigar::Equal(10),
            Cigar::Ins(2),
            Cigar::Equal(6),
            Cigar::Diff(4),
        ];

        let calculator = RangeIdentityCalculator::new(&CigarString(cigar_str));
        println!("{:?}", calculator);

        // 测试完全覆盖的区域
        assert!(
            (calculator.range_identity(0, 5) - 1.0).abs() < 1e-6,
            "{}",
            calculator.range_identity(0, 5)
        ); // 全是 Match
        assert!(
            (calculator.range_identity(5, 8) - 0.0).abs() < 1e-6,
            "{}",
            calculator.range_identity(5, 8)
        ); // 全是 Diff
        assert!(
            (calculator.range_identity(8, 18) - 1.0).abs() < 1e-6,
            "{}",
            calculator.range_identity(8, 18)
        ); // 全是 Equal

        // 测试部分覆盖的区域
        assert!(
            (calculator.range_identity(3, 7) - 0.5).abs() < 1e-6,
            "{}",
            calculator.range_identity(3, 7)
        ); // Match 和 Diff 混合
        assert!(
            (calculator.range_identity(7, 12) - 0.8).abs() < 1e-6,
            "{}",
            calculator.range_identity(7, 12)
        ); // Diff 和 Equal 混合

        // 测试包含 Insert 的区域
        assert!(
            (calculator.range_identity(18, 26) - 0.75).abs() < 1e-6,
            "{}",
            calculator.range_identity(18, 26)
        ); // Equal + Ins + Equal
    }

    #[test]
    fn test_edge_cases() {
        let cigar_str = vec![Cigar::Equal(5), Cigar::Diff(5)];

        let calculator = RangeIdentityCalculator::new(&CigarString(cigar_str));

        // 测试完全在 Match 的边界
        assert!(
            (calculator.range_identity(0, 3) - 1.0).abs() < 1e-6,
            "{}",
            calculator.range_identity(0, 3)
        );

        // 测试完全在 Diff 的边界
        assert!(
            (calculator.range_identity(6, 10) - 0.0).abs() < 1e-6,
            "{}",
            calculator.range_identity(6, 10)
        );

        // 跨越 Match 和 Diff
        assert!(
            (calculator.range_identity(4, 7) - 0.333333).abs() < 1e-6,
            "{}",
            calculator.range_identity(4, 7)
        );
    }


    #[test]
    fn test_parse_cigar_str() {
        println!("{:?}", parse_cigar_string("4=3S"));
    }
}
