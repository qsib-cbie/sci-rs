use core::{borrow::Borrow, cmp::min, iter::Sum, ops::SubAssign};
use nalgebra::{ClosedAdd, ClosedMul, DVector, Scalar};
use num_traits::{Float, One, Zero};

use super::{design::Sos, pad, sosfilt_st, sosfilt_zi_dyn, Pad};

///
/// A forward-backward digital filter using cascaded second-order sections
///
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfiltfilt.html#scipy.signal.sosfiltfilt>
///
///
pub fn sosfiltfilt_dyn<YI, F, const N: usize>(y: YI, sos: &[Sos<F>; N]) -> Vec<F>
where
    F: Float + PartialEq + Scalar + Zero + One + ClosedMul + ClosedAdd + Sum + SubAssign,
    YI: Iterator,
    YI::Item: Borrow<F>,
{
    let ntaps = 2 * N + 1;
    let bzeros = sos
        .iter()
        .map(|s| if s.b[2] == F::zero() { 1 } else { 0 })
        .sum::<usize>();
    let azeros = sos
        .iter()
        .map(|s| if s.a[2] == F::zero() { 1 } else { 0 })
        .sum::<usize>();
    let ntaps = ntaps - min(bzeros, azeros);
    let y = y.map(|yi| *yi.borrow()).collect::<Vec<F>>();
    let y_len = y.len();
    let x = DVector::<F>::from_vec(y);
    let (edge, ext) = pad(Pad::Odd, None, x, 0, ntaps);

    let mut init_sos = sos.clone();
    sosfilt_zi_dyn::<_, _, Sos<F>>(init_sos.iter_mut());

    let x0 = *ext.index(0);
    let mut sos_x = init_sos.clone();
    for s in sos_x.iter_mut() {
        s.zi0 *= x0;
        s.zi1 *= x0;
    }
    let y = sosfilt_st(ext.iter(), &sos_x).collect::<Vec<_>>();

    let y0 = *y.last().unwrap();
    let mut sos_y = init_sos.clone();
    for s in sos_y.iter_mut() {
        s.zi0 *= y0;
        s.zi1 *= y0;
    }
    let mut z = sosfilt_st(y.iter().rev(), &sos_y)
        .skip(edge)
        .take(y_len)
        .collect::<Vec<_>>();
    z.reverse();
    z
}

#[cfg(test)]
mod tests {
    use dasp::{signal, Signal};
    use gnuplot::Figure;
    use itertools::Itertools;

    use super::*;

    #[test]
    fn can_sosfiltfilt() {
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
        // println!("{:?}", &sin_wave);

        let bp_wave = sosfiltfilt_dyn(sin_wave.iter(), &sos);
        println!("{:?}", bp_wave);
        assert_eq!(sin_wave.len(), bp_wave.len());

        let mut fig = Figure::new();

        fig.axes2d().lines(
            bp_wave
                .iter()
                .take(400)
                .enumerate()
                .map(|(i, _)| i)
                .collect_vec(),
            bp_wave.iter().take(400).collect_vec(),
            &[],
        );

        fig.set_pre_commands("set term dumb 100 30");
        fig.show().unwrap();

        println!("{:?}", &bp_wave[..10]);
        println!("{:?}", &sin_wave[..10]);
    }

    #[test]
    fn can_sosfiltfilt_f32() {
        // 4th order butterworth bandpass 10 to 50 at 1666Hz
        let filter: [f32; 24] = [
            2.6775767382597835e-05 as f32,
            5.355153476519567e-05 as f32,
            2.6775767382597835e-05 as f32,
            1.0 as f32,
            -1.7991202154617734 as f32,
            0.8162578614819005 as f32,
            1.0 as f32,
            2.0 as f32,
            1.0 as f32,
            1.0 as f32,
            -1.8774769894419825 as f32,
            0.9094302413068086 as f32,
            1.0 as f32,
            -2.0 as f32,
            1.0 as f32,
            1.0 as f32,
            -1.9237959892866103 as f32,
            0.9263794671616161 as f32,
            1.0 as f32,
            -2.0 as f32,
            1.0 as f32,
            1.0 as f32,
            -1.978497311228862 as f32,
            0.9799894886973378 as f32,
        ];
        let sos = Sos::<f32>::from_scipy::<24, 4>(filter);

        // A signal with a frequency that we can recover
        let sample_hz = 1666.;
        let seconds = 10;
        let mut signal = signal::rate(sample_hz).const_hz(25.).sine();
        let sin_wave: Vec<f32> = (0..seconds * sample_hz as usize)
            .map(|_| signal.next() as f32)
            .collect_vec();
        // println!("{:?}", &sin_wave);

        let bp_wave = sosfiltfilt_dyn(sin_wave.iter(), &sos);
        assert_eq!(sin_wave.len(), bp_wave.len());
        println!("{:?}", bp_wave);

        let mut fig = Figure::new();

        fig.axes2d().lines(
            bp_wave
                .iter()
                .take(400)
                .enumerate()
                .map(|(i, _)| i)
                .collect_vec(),
            bp_wave.iter().take(400).collect_vec(),
            &[],
        );

        fig.set_pre_commands("set term dumb 100 30");
        fig.show().unwrap();

        println!("{:?}", &bp_wave[..10]);
        println!("{:?}", &sin_wave[..10]);
    }
}
