pub fn phreq_list_2_quality(baseq: &[u8]) -> Option<f32> {
    if baseq.len() == 0 {
        return None;
    }
    let length = baseq.len() as f64;
    let err_rate = baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .reduce(|acc, v| acc + v)
        .unwrap()
        / length;
    Some((1.0_f64 - err_rate) as f32)
}

pub fn phreq_list_2_quality_list(baseq: &[u8]) -> Vec<f32> {
    baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .map(|v| 1.0 - v)
        .map(|v| v as f32)
        .collect()
}

pub fn phreq_list_2_error_list(baseq: &[u8]) -> Vec<f32> {
    baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .map(|v| v as f32)
        .collect()
}

pub fn phreq2err(phreq: f64) -> f64 {
    10.0_f64.powf(phreq / -10.0_f64)
}

pub fn phreq2quality(phreq: f64) -> f64 {
    1.0 - phreq2err(phreq)
}

pub fn quality_2_phreq_f32(mut quality: f32, eps: Option<f32>) -> f32 {
    let eps = eps.unwrap_or(1e-5);
    let max_quality = 1.0_f32 - eps;

    quality = if quality > max_quality {
        max_quality
    } else {
        quality
    };
    -10.0_f32 * (1.0_f32 - quality).log10()
}

pub fn quality_2_phreq(quality: f32, eps: Option<f32>) -> u8 {
    let phreq = quality_2_phreq_f32(quality, eps);
    phreq.round() as u8
}

#[cfg(test)]
mod test {
    use crate::phreq::phreq_list_2_error_list;

    use super::{phreq2quality, phreq_list_2_quality, phreq_list_2_quality_list, quality_2_phreq};

    #[test]
    fn test_phreq_list_2_quality() {
        let quality = phreq_list_2_quality(&[20, 20, 20, 20]);
        assert!((quality.unwrap() - 0.99) < 1e-3);

        let quality = phreq_list_2_quality(&[30, 30, 30, 30]);
        assert!((quality.unwrap() - 0.999) < 1e-4);
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

    #[test]
    fn test_phreq_list_2_quality_list() {
        let quality = phreq_list_2_quality_list(&[20, 30, 40, 50]);
        assert!((quality[0] - 0.99).abs() < 1e-6);
        assert!((quality[1] - 0.999).abs() < 1e-6);
        assert!((quality[2] - 0.9999).abs() < 1e-6);
        assert!((quality[3] - 0.99999).abs() < 1e-6);
    }

    #[test]
    fn test_phreq_list_2_err_list() {
        let quality = phreq_list_2_error_list(&[20, 30, 40, 50]);
        assert!((quality[0] - (1.0 - 0.99)).abs() < 1e-6);
        assert!((quality[1] - (1.0 - 0.999)).abs() < 1e-6);
        assert!((quality[2] - (1.0 - 0.9999)).abs() < 1e-6);
        assert!((quality[3] - (1.0 - 0.99999)).abs() < 1e-6);
    }

    #[test]
    fn test_phreq2quality() {
        let quality = phreq2quality(20.0);
        assert!((quality - 0.99).abs() < 1e-6);
    }
}
