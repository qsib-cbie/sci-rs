use core::borrow::Borrow;
use itertools::Itertools;
use num_traits::Float;

// #[cfg(feature = "use_std")]
// use itertools::Itertools;

///
/// Second Order Section Representation
/// of a digital filter transfer function
///
/// Example as Biquad Digital filter from
/// https://en.wikipedia.org/wiki/Digital_biquad_filter
///
/// H(z) =
///        (b0 + b1 * z_1 + b2 * z_2) /
///        (a0 + a1 * z_1 + a2 * z_2)
///
/// Sos { b: [b0, b1, b2], a: [a0, a1, a2] }
///
#[derive(Debug, Copy, Clone)]
pub struct Sos<F: Float> {
    ///
    /// Transfer coefficients
    ///
    b: [F; 3],
    a: [F; 3],

    ///
    /// Filter Delay Values
    ///
    zi0: F,
    zi1: F,
}

impl<F: Float> Default for Sos<F> {
    fn default() -> Self {
        Self {
            b: [F::zero(); 3],
            a: [F::zero(); 3],
            zi0: F::zero(),
            zi1: F::zero(),
        }
    }
}

impl<F: Float> Sos<F> {
    // Normalizes to Direct Form 1 where a[0] is 1
    pub fn new(mut b: [F; 3], mut a: [F; 3]) -> Sos<F> {
        b[0] = b[0] / a[0];
        b[1] = b[1] / a[0];
        b[2] = b[2] / a[0];
        a[0] = a[0] / a[0];
        a[1] = a[1] / a[0];
        a[2] = a[2] / a[0];
        Sos::<F> {
            b,
            a,
            ..Default::default()
        }
    }

    pub fn from_scipy<const M: usize, const N: usize>(sos: [F; M]) -> [Sos<F>; N] {
        // TODO: Replace with stable const expressions
        assert!(N * 6 == M);

        let mut rslt = [Default::default(); N];
        sos.iter()
            .tuples::<(&F, &F, &F, &F, &F, &F)>()
            .take(N)
            .zip(rslt.iter_mut())
            .for_each(|(ba, sos)| {
                *sos = Sos::new([*ba.0, *ba.1, *ba.2], [*ba.3, *ba.4, *ba.5]);
            });
        rslt
    }
}

///
/// A series of Second Order Sections may be used to
/// filter a stream of inputs with a normalization factor
///
pub struct SosFilt<I, F: Float, const N: usize> {
    pub(crate) iter: I,

    ///
    /// H(z) = product(H(z) for each Sos)
    ///
    sos: [Sos<F>; N],

    /// Normalization factor on the output H(z)
    ///
    g: F,
}

impl<I, B, F, const N: usize> Iterator for SosFilt<I, F, N>
where
    I: Iterator<Item = B>,
    B: Borrow<F>,
    F: Float,
{
    type Item = F;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|xn| {
            // Based on sosfilt: https://github.com/scipy/scipy/blob/v1.7.1/scipy/signal/_sosfilt.pyx
            self.sos.iter_mut().fold(*xn.borrow(), |x_cur, sos| {
                let x_new = sos.b[0] * x_cur + sos.zi0;
                sos.zi0 = sos.b[1] * x_cur - sos.a[1] * x_new + sos.zi1;
                sos.zi1 = sos.b[2] * x_cur - sos.a[2] * x_new;
                x_new
            }) * self.g
        })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

///
/// Second Order Sections filter an iterator
///
pub fn sosfilt<YI, F, const N: usize>(y: YI, sos: &[Sos<F>; N]) -> SosFilt<YI, F, N>
where
    F: Float,
    YI: Iterator,
    YI::Item: Borrow<F>,
{
    SosFilt {
        iter: y,
        sos: sos.clone(),
        g: F::from(1.0).unwrap(),
    }
}

#[cfg(test)]
mod tests {
    use dasp::{signal, Signal};
    use gnuplot::Figure;

    use super::*;

    #[test]
    fn can_use_external_filter_design() {
        let _design_filter = r#"
import scipy.signal as sg
buttersos = sg.butter(4, [0.5, 100], btype='bandpass', output='sos', fs=1666)
print(buttersos)
        "#;

        let butterworth_filter_sos = [
            0.00474269,
            0.00948539,
            0.00474269,
            1.,
            -1.05531479,
            0.29986557,
            1.,
            2.,
            1.,
            1.,
            -1.32397785,
            0.6355536,
            1.,
            -2.,
            1.,
            1.,
            -1.99416225,
            0.99417226,
            1.,
            -2.,
            1.,
            1.,
            -1.99760501,
            0.9976149,
        ];
        let sos: [Sos<f64>; 4] = Sos::from_scipy(butterworth_filter_sos);
        println!("{:?}", sos);
    }

    #[test]
    fn can_iter() {
        let mut count = 0;
        let sos = [Sos::<f64> {
            b: [0.00474269, 0.00948539, 0.00474269],
            a: [1., -1.05531479, 0.29986557],
            ..Default::default()
        }];

        const NUM_ITERS: usize = 100;
        sosfilt((0..NUM_ITERS).map(|i| i as f64), &sos).for_each(|_| {
            count += 1;
        });

        assert_eq!(NUM_ITERS, count);
    }

    #[test]
    fn can_bandpass() {
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
        println!("{:?}", &sin_wave);

        let bp_wave = sosfilt(sin_wave.iter(), &sos).collect_vec();
        // println!("{:?}", bp_wave);

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
