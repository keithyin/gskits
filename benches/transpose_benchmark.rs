use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gskits::matrix::transpose;

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

criterion_group!(benches, transpose_benchmark, transpose_blocked_benchmark);
criterion_main!(benches);
