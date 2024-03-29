use criterion::{black_box, criterion_group, criterion_main, Criterion};
use dasp_signal::{rate, Signal};
use sci_rs::signal::filter::design::Sos;
use sci_rs::signal::filter::{sosfilt_dyn, sosfilt_fast32_st};

// TLDR: 8.5x faster
// sosfilt_st is as fast as sosfilt_dyn

///
/// 4th order Butterworth Bandpass Sosfilt 10 seconds of 1666Hz sine wave
///
/// Python scipy.signal.sosfilt:
/// ```
/// avg over 1000 runs 66,776,329 ns/iter
/// ```
///
/// Rust implementation
/// ```
/// sosfilt_100x_dyn        time:   [7.7290 ms 7.7425 ms 7.7572 ms]
/// ```
///
fn butter_sosfilt_100x_dyn(c: &mut Criterion) {
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
    let mut sos = Sos::from_scipy_dyn(4, filter.to_vec());

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next())
        .collect::<Vec<_>>();
    let sin_wave = (0..100).flat_map(|_| sin_wave.clone()).collect::<Vec<_>>();

    c.bench_function("sosfilt_100x_dyn", |b| {
        b.iter(|| {
            black_box(sosfilt_dyn(sin_wave.iter(), &mut sos));
        });
    });
}

fn butter_sosfilt_f64(c: &mut Criterion) {
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
    let mut sos = Sos::from_scipy_dyn(4, filter.to_vec());

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next())
        .collect::<Vec<_>>();

    c.bench_function("sosfilt_f64", |b| {
        b.iter(|| {
            black_box(sosfilt_dyn(sin_wave.iter(), &mut sos));
        });
    });
}

fn butter_sosfilt_f32(c: &mut Criterion) {
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
    let mut sos = Sos::from_scipy_dyn(4, filter.to_vec());

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f32> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next() as f32)
        .collect::<Vec<_>>();

    c.bench_function("sosfilt_f32", |b| {
        b.iter(|| {
            black_box(sosfilt_dyn(sin_wave.iter(), &mut sos));
        });
    });
}

///
/// 4th order Butterworth Bandpass Sosfilt 10 seconds of 1666Hz sine wave
/// 2x faster than sosfilt_dyn due to tiling and cpu pipelining
///
/// ```
/// sosfilt_f32             time:   [139.61 µs 139.86 µs 140.16 µs]
/// sosfilt_fast32_st4       time:   [78.412 µs 79.758 µs 81.399 µs]
/// ```
///
///
fn butter_sosfilt_fast32_st4(c: &mut Criterion) {
    // 4th order butterworth bandpass 10 to 50 at 1666Hz
    let filter: [f32; 24] = [
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
    let mut sos = Sos::from_scipy_dyn(4, filter.to_vec());

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f32> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next() as f32)
        .collect::<Vec<_>>();
    let mut buf = vec![0.0; sin_wave.len()];

    c.bench_function("sosfilt_fast32_st4", |b| {
        b.iter(|| {
            black_box(sosfilt_fast32_st(&sin_wave, &mut sos, &mut buf));
        });
    });
}

///
/// 4th vs 8th order Butterworth Bandpass Sosfilt 10 seconds of 1666Hz sine wave
/// 2x faster due to tiling and cpu pipelining
///
/// ```
/// // with tiling specialization
/// sosfilt_fast32_st4       time:   [78.412 µs 79.758 µs 81.399 µs]
/// sosfilt_fast32_st8       time:   [133.39 µs 133.73 µs 134.08 µs]
///
/// // without tiling specialization
/// sosfilt_fast32_st8      time:   [536.05 µs 537.17 µs 538.31 µs]
///    change: [+300.64% +301.99% +303.35%] (p = 0.00 < 0.05)
///    Performance has regressed.
/// ```
///
///
fn butter_sosfilt_fast32_st8(c: &mut Criterion) {
    // 8th order butterworth bandpass 10 to 50 at 1666HzA
    let filter: [f32; 48] = [
        7.223657016655901e-10,
        1.4447314033311803e-09,
        7.223657016655901e-10,
        1.0,
        -1.8117367715812775,
        0.82390242865626,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.8027320089797896,
        0.8243217191895463,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.8436432946399057,
        0.8727119112095707,
        1.0,
        2.0,
        1.0,
        1.0,
        -1.898284650388357,
        0.9018968456559892,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.9415418844285526,
        0.943617934088245,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.918327282768592,
        0.9524499384015938,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.967653721895735,
        0.9692546426434141,
        1.0,
        -2.0,
        1.0,
        1.0,
        -1.9886761584321624,
        0.9901117398066808,
    ];
    let mut sos = Sos::from_scipy_dyn(8, filter.to_vec());

    // A signal with a frequency that we can recover
    let sample_hz = 1666.;
    let seconds = 10;
    let mut signal = rate(sample_hz).const_hz(25.).sine();
    let sin_wave: Vec<f32> = (0..seconds * sample_hz as usize)
        .map(|_| signal.next() as f32)
        .collect::<Vec<_>>();
    let mut buf = vec![0.0; sin_wave.len()];

    c.bench_function("sosfilt_fast32_st8", |b| {
        b.iter(|| {
            black_box(sosfilt_fast32_st(&sin_wave, &mut sos, &mut buf));
        });
    });
}

criterion_group!(
    benches,
    butter_sosfilt_100x_dyn,
    butter_sosfilt_f64,
    butter_sosfilt_f32,
    butter_sosfilt_fast32_st4,
    butter_sosfilt_fast32_st8,
);
criterion_main!(benches);
