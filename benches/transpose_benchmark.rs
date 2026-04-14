use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gskits::matrix::transpose::{self, avx2_transopose_u8_16x16};

fn transpose_benchmark(c: &mut Criterion) {
    let values = (0..100000).into_iter().collect::<Vec<_>>();
    c.bench_function("transpose", |b| {
        b.iter(|| transpose::transpose(black_box(&values), black_box(10), black_box(10000)))
    });
}

fn transpose_blocked_benchmark(c: &mut Criterion) {
    let values = (0..(3375_usize * 8870))
        .into_iter()
        .map(|v| (v % 256) as u8)
        .collect();

    let mut group = c.benchmark_group("transpose_blocked_benchmark");
    for block_size in [
        None,
        Some(8),
        Some(16),
        Some(32),
        Some(64),
        Some(128),
        Some(256),
        Some(512),
    ] {
        group.bench_function(format!("{:?}", block_size), |b| {
            b.iter(|| {
                transpose::transpose_blocked(
                    black_box(&values),
                    black_box(3375),
                    black_box(8870),
                    black_box(block_size),
                )
            })
        });
    }
}

fn transpose_benchmark_transpose_16x16(c: &mut Criterion) {
    let value = (0..(16 * 16))
        .into_iter()
        .map(|v| v as u8)
        .collect::<Vec<u8>>();

    let mut output = vec![0; 16 * 16];

    c.bench_function("transpose_16x16", |b| {
        b.iter(|| {
            transpose::transpose_withoutput(
                black_box(&value),
                black_box(&mut output),
                black_box(16),
                black_box(16),
            );
        })
    });
}

fn transpose_benchmark_transpose_16x16_avx2(c: &mut Criterion) {
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

    c.bench_function("transpose_16x16_avx2", |b| {
        b.iter(|| {
            transpose::avx2_transopose_u8_16x16(black_box(&value2), black_box(&mut output2));
        })
    });
}

criterion_group!(
    benches,
    transpose_benchmark,
    transpose_blocked_benchmark,
    transpose_benchmark_transpose_16x16,
    transpose_benchmark_transpose_16x16_avx2
);
criterion_main!(benches);
