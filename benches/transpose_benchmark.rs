use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gskits::matrix::transpose;

fn transpose_benchmark(c: &mut Criterion) {
    let values = (0..100000).into_iter().collect::<Vec<_>>();
    c.bench_function("transpose", |b| {
        b.iter(|| transpose::transpose(black_box(&values), black_box(10), black_box(10000)))
    });
}

fn transpose_blocked_benchmark(c: &mut Criterion) {
    let values = (0..100000).into_iter().collect::<Vec<_>>();
    c.bench_function("transpose_block", |b| {
        b.iter(|| {
            transpose::transpose_blocked(
                black_box(&values),
                black_box(10),
                black_box(10000),
                black_box(None),
            )
        })
    });
    
}

criterion_group!(benches, transpose_benchmark, transpose_blocked_benchmark);
criterion_main!(benches);
