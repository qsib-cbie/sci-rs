use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::{array, Array1, ArrayView1};
use rand::rngs::ThreadRng;
use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
use sci_rs::signal::filter::LFilter;
use sci_rs::signal::windows::Hamming;
use std::num::NonZeroUsize;

/// Get a randomized signal from instance of `rng`.
fn randomized_signal(
    mut rng: ThreadRng,
    num_freqs: NonZeroUsize,
    num_data_points: NonZeroUsize,
    time_seconds: f64,
) -> (Array1<f64>, Array1<f64>) {
    use rand::RngExt;

    let nf: usize = num_freqs.into(); // Num freqs
    let n: usize = num_data_points.into(); // Num data points
    let num_secs: f64 = time_seconds; // this corresponds to how much N data points corresponds to in time.
    let nyq_freq: f64 = 0.5 * n as f64 / num_secs; // Do not generated frequency should be larger than this

    let t: Array1<f64> = Array1::linspace(0.0, num_secs, n); // Time-axis
    let global_shift = t.mapv(|ti| 7.0 * (-ti / 2.0).exp()); // Total offset independent of sinusoidal generation

    let ampl = {
        let orig: Array1<f64> = Array1::from_iter((0..nf).map(|_| rng.random_range(0.5..1.5)));
        let decay: Array1<f64> = Array1::from_iter((0..nf).map(|i| 1. / (1.1f64.powf(i as f64)))); // Weight

        decay * orig
    };
    let freqs = {
        let mut freqs = Vec::with_capacity(nf);
        const INITIAL_FREQ: f64 = 1.2;
        freqs.push(INITIAL_FREQ);
        (1..nf).fold(INITIAL_FREQ, |acc, _| {
            let next_freq = (acc * rng.random_range(1.01..2.01)) % nyq_freq; // Approximately double wrt the previous frequency.
            freqs.push(next_freq);
            next_freq
        });
        freqs
    };
    let phases: Vec<_> = (0..nf)
        .map(|_| rng.random_range(0.0..std::f64::consts::PI))
        .collect();

    let mut result: Array1<f64> = Array1::zeros((n,)) + global_shift;

    for ((a, freq), p) in ampl
        .into_iter()
        .zip(freqs.into_iter())
        .zip(phases.into_iter())
    {
        let wave = t.mapv(|ti| a * (freq * ti + p).sin());
        result += &wave;
    }

    (t, result)
}

/// Test the filter with zi.
///
/// Use window asasociated with `decimate`'s default values, running at decimation factor = 500.
fn lfilter_dyn(c: &mut Criterion) {
    const DECIMATION_FACTOR: usize = 50;
    const FILTER_ORDER: usize = DECIMATION_FACTOR * 20;

    // Finite impulse response from Hamming window;
    let b: Array1<f64> = firwin_dyn(
        FILTER_ORDER + 1,
        &[1. / (DECIMATION_FACTOR as f64)],
        None,
        None::<&Hamming>,
        &FilterBandType::Lowpass,
        None,
        None,
    )
    .unwrap()
    .into();
    let a = array![1.];

    let (_, signal) = randomized_signal(
        rand::rng(),
        NonZeroUsize::new(14).unwrap(),
        NonZeroUsize::new(1 << 16).unwrap(),
        15.,
    );

    // Apply with criterion
    c.bench_with_input(
        BenchmarkId::new("lfilter_dyn", DECIMATION_FACTOR),
        &signal,
        |bench, sig| {
            bench.iter(|| {
                ArrayView1::lfilter(
                    black_box((&b).into()),
                    black_box((&a).into()),
                    black_box((sig).into()),
                    None,
                    None,
                )
            })
        },
    );
}

criterion_group!(benches, lfilter_dyn);
criterion_main!(benches);
