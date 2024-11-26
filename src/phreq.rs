pub fn phreq_list_2_quality(baseq: &[u8]) -> f32 {
    let length = baseq.len() as f64;
    let err_rate = baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .reduce(|acc, v| acc + v)
        .unwrap()
        / length;
    (1.0_f64 - err_rate) as f32
}

pub fn quality_2_phreq(mut quality: f32, eps: Option<f32>) -> u8 {
    let eps = eps.unwrap_or(1e-5);
    let max_quality = 1.0_f32 - eps;

    quality = if quality > max_quality {
        max_quality
    } else {
        quality
    };
    let phreq = -10.0_f32 * (1.0_f32 - quality).log10();
    phreq.round() as u8
}


#[cfg(test)]
mod test {
    use super::{phreq_list_2_quality, quality_2_phreq};


    #[test]
    fn test_phreq_list_2_quality() {
        let quality = phreq_list_2_quality(&[20, 20, 20, 20]);
        assert!((quality - 0.99) < 1e-3);

        let quality = phreq_list_2_quality(&[30, 30, 30, 30]);
        assert!((quality - 0.999) < 1e-4);
    }

    #[test]
    fn test_quality_2_phreq() {

        let phreq = quality_2_phreq(0.99, None);
        assert_eq!(phreq, 20);

        let phreq = quality_2_phreq(0.999, None);
        assert_eq!(phreq, 30);

        let phreq = quality_2_phreq(0.9999, None);
        assert_eq!(phreq, 40);

    }
}