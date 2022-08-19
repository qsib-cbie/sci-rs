use core::borrow::Borrow;
use num_traits::Float;

use super::design::Sos;

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
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfilt.html>
///
/// `sos` holds the `zf` return value from scipy, so reusing a `sos` reference
/// from a previous iteration achieves the same result as passing `zi` in the scipy interface
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
    use itertools::Itertools;

    use super::*;

    #[test]
    fn can_sosfilt() {
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
