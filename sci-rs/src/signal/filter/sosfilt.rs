use core::{
    borrow::Borrow,
    mem::{transmute, MaybeUninit},
};
use nalgebra::RealField;
use num_traits::Float;

use super::design::Sos;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

///
/// A series of Second Order Sections may be used to
/// filter a stream of inputs with a normalization factor
///
#[cfg(feature = "alloc")]
#[inline]
///
/// Second Order Sections filter an iterator
///
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfilt.html>
///
/// `sos` holds the `zf` return value from scipy, so reusing a `sos` reference
/// from a previous iteration achieves the same result as passing `zi` in the scipy interface
///
pub fn sosfilt_dyn<YI, F>(y: YI, sos: &mut [Sos<F>]) -> Vec<F>
where
    F: RealField + Copy,
    YI: IntoIterator,
    YI::Item: Borrow<F>,
{
    match sos.len() {
        2 => y
            .into_iter()
            .map(|yi0| {
                sos[..2].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        3 => y
            .into_iter()
            .map(|yi0| {
                sos[..3].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        4 => y
            .into_iter()
            .map(|yi0| {
                sos[..4].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        5 => y
            .into_iter()
            .map(|yi0| {
                sos[..5].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        6 => y
            .into_iter()
            .map(|yi0| {
                sos[..6].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        8 => y
            .into_iter()
            .map(|yi0| {
                sos[..8].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        10 => y
            .into_iter()
            .map(|yi0| {
                sos[..10].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        12 => y
            .into_iter()
            .map(|yi0| {
                sos[..12].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        14 => y
            .into_iter()
            .map(|yi0| {
                sos[..14].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        16 => y
            .into_iter()
            .map(|yi0| {
                sos[..16].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        18 => y
            .into_iter()
            .map(|yi0| {
                sos[..18].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        20 => y
            .into_iter()
            .map(|yi0| {
                sos[..20].iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                    let x_new = sos.b[0] * x_curr + sos.zi0;
                    sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                    sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                    x_new
                })
            })
            .collect::<Vec<_>>(),
        _ => y
            .into_iter()
            .map(|yi0| sosfilt_item(*yi0.borrow(), sos))
            .collect::<Vec<_>>(),
    }
}

///
/// Apply the cascaded Biquad filter represented by `sos` to the input `y`
/// representing a single sample. This avoids allocating at the cost of not
/// being able to specialize the filter length.
///
/// This is nearly always slower for long iterators than `sosfilt_dyn` due to
/// the overhead of the `sos` slice iteration and lack of pipelining.
///
pub fn sosfilt_st<'it, YI, F>(y: YI, sos: &'it mut [Sos<F>]) -> impl Iterator<Item = F> + 'it
where
    F: RealField + Copy,
    YI: IntoIterator + 'it,
    YI::Item: Borrow<F>,
{
    y.into_iter().map(|yi0| sosfilt_item(*yi0.borrow(), sos))
}

///
/// Apply the cascaded Biquad filter represented by `sos` to the input `y`
/// representing a single sample. This may be useful if the filter is being
/// applied to a single sample at a time as data is discretely sampled.
///
/// This is nearly always slower than `sosfilt_dyn` due to the overhead of
/// the `sos` slice iteration and lack of pipelining.
///
#[inline(always)]
pub fn sosfilt_item<F, B>(y: B, sos: &mut [Sos<F>]) -> F
where
    F: RealField + Copy,
    B: Borrow<F>,
{
    sos.iter_mut().fold(*y.borrow(), |x_curr, sos| {
        let x_new = sos.b[0] * x_curr + sos.zi0;
        sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
        sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
        x_new
    })
}

type Sos32 = Sos<f32>;

#[inline(always)]
fn biquad_fold(yi: f32, sos: &mut Sos32) -> f32 {
    let x = sos.b[0] * yi + sos.zi0;
    sos.zi0 = sos.b[1] * yi - sos.a[1] * x + sos.zi1;
    sos.zi1 = sos.b[2] * yi - sos.a[2] * x;
    x
}

fn _sosfilt32(y: &[f32], sos: &mut [Sos32], z: &mut [f32]) {
    if y.len() != z.len() {
        panic!();
    }
    if y.is_empty() {
        return;
    }
    const TILE: usize = 4;
    match (sos.len(), y.len() % TILE == 0) {
        (2, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = *yi;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt32(&y[idx..], sos, &mut z[idx..]);
        }
        (4, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = *yi;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt32(&y[idx..], sos, &mut z[idx..]);
        }
        (6, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = *yi;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt32(&y[idx..], sos, &mut z[idx..]);
        }
        (8, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = *yi;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt32(&y[idx..], sos, &mut z[idx..]);
        }
        (10, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = *yi;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt32(&y[idx..], sos, &mut z[idx..]);
        }
        _ => {
            for (yi, zi) in y.iter().zip(z.iter_mut()) {
                *zi = *yi;
                for s in sos.iter_mut() {
                    *zi = biquad_fold(*zi, s);
                }
            }
        }
    }
}

fn _sosfilt_isize_32<I: Copy + Into<isize>>(y: &[I], sos: &mut [Sos32], z: &mut [f32]) {
    if y.len() != z.len() {
        panic!();
    }
    if y.is_empty() {
        return;
    }
    const TILE: usize = 4;
    match (sos.len(), y.len() % TILE == 0) {
        (2, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = Into::<isize>::into(*yi) as f32;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt_isize_32(&y[idx..], sos, &mut z[idx..]);
        }
        (4, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = Into::<isize>::into(*yi) as f32;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt_isize_32(&y[idx..], sos, &mut z[idx..]);
        }
        (6, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = Into::<isize>::into(*yi) as f32;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt_isize_32(&y[idx..], sos, &mut z[idx..]);
        }
        (8, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = Into::<isize>::into(*yi) as f32;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt_isize_32(&y[idx..], sos, &mut z[idx..]);
        }
        (10, true) => {
            let rem = y.len() % TILE;
            y.chunks_exact(TILE)
                .zip(z.chunks_exact_mut(TILE))
                .for_each(|c| {
                    for (yi, zi) in c.0.iter().zip(c.1.iter_mut()) {
                        *zi = Into::<isize>::into(*yi) as f32;
                        for s in sos.iter_mut() {
                            *zi = biquad_fold(*zi, s);
                        }
                    }
                });
            let idx = y.len() - rem;
            _sosfilt_isize_32(&y[idx..], sos, &mut z[idx..]);
        }
        _ => {
            for (yi, zi) in y.iter().zip(z.iter_mut()) {
                *zi = Into::<isize>::into(*yi) as f32;
                for s in sos.iter_mut() {
                    *zi = biquad_fold(*zi, s);
                }
            }
        }
    }
}

///
/// A specialized cascaded Biquad filter for 32-bit floating point samples
///
/// Including acceleratored
///  * Single-sided 4th and 8th order filters
///     * Example: 4th or 8th order lowpass Butterworth
///  * Double-sided 4th order filters are accelerated
///     * Example: 4th order bandpass Butterworth
///
pub fn sosfilt_fast32_st(y: &[f32], sos: &mut [Sos32], z: &mut [f32]) {
    _sosfilt32(y, sos, z);
}

///
/// A specialized cascaded Biquad filter for signed samples
/// filtered as 32-bit floating point samples. Samples implement
/// Into<isize> to convert to a signed integer for speed on native bitwidths.
///
/// Including acceleratored
///  * Single-sided 4th and 8th order filters
///     * Example: 4th or 8th order lowpass Butterworth
///  * Double-sided 4th order filters are accelerated
///     * Example: 4th order bandpass Butterworth
///
pub fn sosfilt_ifast32_st<I>(y: &[I], sos: &mut [Sos32], z: &mut [f32])
where
    I: Into<isize> + Copy,
{
    _sosfilt_isize_32(y, sos, z);
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use dasp_signal::{rate, Signal};

    #[cfg(feature = "alloc")]
    #[test]
    fn can_sosfilt() {
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
        assert_eq!(sos.len(), 4);

        // A signal with a frequency that we can recover
        let sample_hz = 1666.;
        let seconds = 10;
        let mut signal = rate(sample_hz).const_hz(25.).sine();
        let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
            .map(|_| signal.next())
            .collect::<Vec<_>>();
        println!("{:?}", &sin_wave);

        let mut sos_item = sos.clone();
        let mut bp_item_wave = Vec::new();
        for yi0 in sin_wave.iter() {
            bp_item_wave.push(sosfilt_item(*yi0, &mut sos_item));
        }
        let mut sos_st = sos.clone();
        let bp_wave = sosfilt_st(sin_wave.iter(), &mut sos_st).collect::<Vec<_>>();
        let mut sos_dyn = sos;
        let bp_dyn_wave = sosfilt_dyn(sin_wave.iter(), &mut sos_dyn);
        for (a, b) in bp_item_wave.iter().zip(bp_wave.iter()) {
            assert_relative_eq!(*a, *b);
        }
        for (a, b) in bp_wave.iter().zip(bp_dyn_wave.iter()) {
            assert_relative_eq!(*a, *b);
        }
        // println!("{:?}", bp_wave);

        println!("{:?}", &bp_wave[..10]);
        println!("{:?}", &sin_wave[..10]);
    }

    #[cfg(feature = "alloc")]
    #[test]
    fn can_resume_sosfilt() {
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
        assert_eq!(sos.len(), 4);

        // A signal with a frequency that we can recover
        let sample_hz = 1666.;
        let seconds = 10;
        let mut signal = rate(sample_hz).const_hz(25.).sine();
        let sin_wave: Vec<f64> = (0..seconds * sample_hz as usize)
            .map(|_| signal.next())
            .collect::<Vec<_>>();
        println!("{:?}", &sin_wave);

        let mut sos_st = sos.clone();
        let bp_wave = sosfilt_dyn(sin_wave.iter(), &mut sos_st);
        let mut sos_st2 = sos;
        let mut bp_st2_wave = sosfilt_dyn(sin_wave.iter().take(sin_wave.len() / 2), &mut sos_st2);
        bp_st2_wave.extend(sosfilt_dyn(
            sin_wave.iter().skip(sin_wave.len() / 2),
            &mut sos_st2,
        ));
        for (a, b) in bp_wave.iter().zip(bp_st2_wave.iter()) {
            assert_relative_eq!(*a, *b);
        }
    }
}
