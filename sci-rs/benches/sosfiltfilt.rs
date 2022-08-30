use criterion::{black_box, criterion_group, criterion_main, Criterion};
use dasp_signal::{rate, Signal};
use sci_rs::signal::filter::{design::Sos, sosfiltfilt_dyn};

/// TLDR: 4.6x faster

///
/// 4th order Butterworth Bandpass Sosfilt 10 seconds of 1666Hz sine wave
///
/// Python scipy.signal.sosfiltfilt:
/// ```
/// avg over 1000 runs 89,924,038 ns/iter
/// ```
///
/// Rust implementation
/// ```
/// sosfiltfilt_100x        time:   [19.412 ms 19.490 ms 19.573 ms]
/// ```
///
fn butter_sosfiltfilt_100x(c: &mut Criterion) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f64; 24] = [
        2.677_576_738_259_783_5e-5,
        5.355_153_476_519_567e-5,
        2.677_576_738_259_783_5e-5,
        1.0,
        -1.7991202154617734,
        0.8162578614819005,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.8774769894419825,
        0.9094302413068086,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.9237959892866103,
        0.9263794671616161,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.978497311228862,
        0.9799894886973378,
    ];
    let sos = Sos::from_scipy::<24, 4>(filter);

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next())
        .collect::<Vec<_>>();
    let sin_wave = (0..100).flat_map(|_| sin_wave.clone()).collect::<Vec<_>>();

    c.bench_function("sosfiltfilt_100x", |b| {
        b.iter(|| {
            black_box(sosfiltfilt_dyn(sin_wave.iter(), &sos));
        });
    });
}

fn butter_sosfiltfilt_10x(c: &mut Criterion) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f64; 24] = [
        2.677_576_738_259_783_5e-5,
        5.355_153_476_519_567e-5,
        2.677_576_738_259_783_5e-5,
        1.0,
        -1.7991202154617734,
        0.8162578614819005,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.8774769894419825,
        0.9094302413068086,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.9237959892866103,
        0.9263794671616161,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.978497311228862,
        0.9799894886973378,
    ];
    let sos = Sos::from_scipy::<24, 4>(filter);

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next())
        .collect::<Vec<_>>();
    let sin_wave = (0..10).flat_map(|_| sin_wave.clone()).collect::<Vec<_>>();

    c.bench_function("sosfiltfilt_10x", |b| {
        b.iter(|| {
            black_box(sosfiltfilt_dyn(sin_wave.iter(), &sos));
        });
    });
}

fn butter_sosfiltfilt_f64(c: &mut Criterion) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f64; 24] = [
        2.677_576_738_259_783_5e-5,
        5.355_153_476_519_567e-5,
        2.677_576_738_259_783_5e-5,
        1.0,
        -1.7991202154617734,
        0.8162578614819005,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.8774769894419825,
        0.9094302413068086,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.9237959892866103,
        0.9263794671616161,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.978497311228862,
        0.9799894886973378,
    ];
    let sos = Sos::from_scipy::<24, 4>(filter);

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next())
        .collect::<Vec<_>>();

    c.bench_function("sosfiltfilt_f64", |b| {
        b.iter(|| {
            black_box(sosfiltfilt_dyn(sin_wave.iter(), &sos));
        });
    });
}

fn butter_sosfiltfilt_f32(c: &mut Criterion) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f32; 24] = [
        2.677_576_8e-5,
        5.355_153_6e-5,
        2.677_576_8e-5,
        1.0,
        -1.799_120_2,
        0.816_257_83,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.877_476_9,
        0.909_430_27,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.923_795_9,
        0.926_379_44,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.978_497_3,
        0.979_989_47,
    ];
    let sos = Sos::from_scipy::<24, 4>(filter);

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f32> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next() as f32)
        .collect::<Vec<_>>();

    c.bench_function("sosfiltfilt_f32", |b| {
        b.iter(|| {
            black_box(sosfiltfilt_dyn(sin_wave.iter(), &sos));
        });
    });
}

criterion_group!(
    benches,
    butter_sosfiltfilt_10x,
    butter_sosfiltfilt_100x,
    butter_sosfiltfilt_f64,
    butter_sosfiltfilt_f32
);
criterion_main!(benches);
