#[macro_use]
extern crate bencher;

use bencher::Bencher;
use dasp::{signal, Signal};
use sci_rs::filt::{sosfilt, Sos};
use itertools::Itertools;

///
/// 4th order Butterworth Bandpass Sosfilt 10 seconds of 1666Hz sine wave
///
/// Python scipy.signal.sosfilt:
/// ```
/// avg over 1000 runs 667.6116660000001 us
/// ```
///
/// Rust implementation
/// ```
/// test butter_sosfilt ... bench:      76,694 ns/iter (+/- 1,336)
/// ```
///
fn butter_sosfilt_f64(bench: &mut Bencher) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f64; 24] = [
        2.6775767382597835e-05,
        5.355153476519567e-05,
        2.6775767382597835e-05,
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
    let mut signal = signal::rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next())
        .collect_vec();

    bench.iter(|| {
        sosfilt(sin_wave.iter(), &sos).collect_vec();
    });
}

fn butter_sosfilt_f32(bench: &mut Bencher) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f32; 24] = [
        2.6775767382597835e-05,
        5.355153476519567e-05,
        2.6775767382597835e-05,
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
    let mut signal = signal::rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f32> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next() as f32)
        .collect_vec();

    bench.iter(|| {
        sosfilt(sin_wave.iter(), &sos).collect_vec();
    });
}

benchmark_group!(benches, butter_sosfilt_f64, butter_sosfilt_f32);
benchmark_main!(benches);
