use core::{
    borrow::Borrow,
    mem::{transmute, MaybeUninit},
};
use nalgebra::RealField;
use num_traits::Float;

use super::design::Sos;
///
/// A series of Second Order Sections may be used to
/// filter a stream of inputs with a normalization factor
///
#[cfg(feature = "use_std")]
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
    YI: Iterator,
    YI::Item: Borrow<F>,
{
    let y = y.collect::<Vec<_>>();
    y.into_iter()
        .map(|yi0| {
            sos.iter_mut().fold(*yi0.borrow(), |x_curr, sos| {
                let x_new = sos.b[0] * x_curr + sos.zi0;
                sos.zi0 = sos.b[1] * x_curr - sos.a[1] * x_new + sos.zi1;
                sos.zi1 = sos.b[2] * x_curr - sos.a[2] * x_new;
                x_new
            })
        })
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use dasp_signal::{rate, Signal};
    #[cfg(feature = "plot")]
    use gnuplot::Figure;

    use super::*;

    #[cfg(feature = "use_std")]
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

        let mut sos_st = sos.clone();
        let bp_wave = sosfilt_dyn(sin_wave.iter(), &mut sos_st);
        let mut sos_dyn = sos;
        let bp_dyn_wave = sosfilt_dyn(sin_wave.iter(), &mut sos_dyn);
        for (a, b) in bp_wave.iter().zip(bp_dyn_wave.iter()) {
            assert_relative_eq!(*a, *b);
        }
        // println!("{:?}", bp_wave);

        #[cfg(feature = "plot")]
        {
            let mut fig = Figure::new();
            fig.axes2d().lines(
                bp_wave
                    .iter()
                    .take(400)
                    .enumerate()
                    .map(|(i, _)| i)
                    .collect::<Vec<_>>(),
                bp_wave.iter().take(400).collect::<Vec<_>>(),
                &[],
            );

            fig.set_pre_commands("set term dumb 100 30");
            fig.show().unwrap();
        }

        println!("{:?}", &bp_wave[..10]);
        println!("{:?}", &sin_wave[..10]);
    }

    #[cfg(feature = "use_std")]
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
