#![no_std]
#[cfg(test)]
#[macro_use]
extern crate std;

use core::{borrow::Borrow, iter::Sum};

use num_traits::Float;

pub fn default<F>() -> F
where
    F: Float + Default,
{
    Default::default()
}

///
/// Autocorrelate the signal `y` with lag `k`,
/// using 1/N formulation
///
/// https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm
///
/// ```
/// use approx::relative_eq;
/// use emb_dsp::autocorr;
///
/// // Flat signal perfectly correlates with itself
/// let y: [f64; 4] = [1.,1.,1.,1.];
/// relative_eq!(1f64, autocorr(y.iter(), 1));
/// ```
///
pub fn autocorr<YI, F>(y: YI, k: usize) -> F
where
    F: Float + Sum,
    YI: Iterator + Clone,
    YI::Item: Borrow<F>,
{
    let y = y.map(|f| f.borrow().clone());
    let n = F::from(y.clone().count()).unwrap();
    let sum: F = y.clone().sum();
    let mean = sum / n;
    let variance: F = y
        .clone()
        .map(|f| {
            let delta = f - mean;
            delta * delta
        })
        .sum::<F>()
        / n;
    let autocovariance: F = y
        .clone()
        .zip(y.skip(k))
        .map(|(fi, fik)| (fi - mean) * (fik - mean))
        .sum::<F>()
        / n;
    autocovariance / variance
}

///
/// Produce an iterator yielding the lag difference, yi1 - yi0,
///
/// https://www.itl.nist.gov/div898/handbook/eda/section3/lagplot.htm
///
/// ```
/// use approx::relative_eq;
/// use emb_dsp::lag_diff;
///
/// // Flat signal perfectly correlates with itself
/// let y: [f64; 4] = [1.,2.,4.,7.];
/// let z = lag_diff(y.iter()).collect::<Vec<_>>();
/// for i in 0..3 {
///     relative_eq!(i as f64 + 1f64, z[i]);
/// }
/// ```
///
pub fn lag_diff<'a, YI, F>(y: YI) -> impl Iterator<Item = F>
where
    F: Float + 'a,
    YI: Iterator + Clone,
    YI::Item: Borrow<F>,
{
    y.clone()
        .zip(y.skip(1))
        .map(|(yi0, yi1)| *yi1.borrow() - *yi0.borrow())
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use std::vec::Vec;

    use approx::assert_relative_eq;
    use gnuplot::{Figure, PlotOption::Caption};
    use itertools::Itertools;

    use super::*;

    #[test]
    fn can_default_float() {
        let rslt = default::<f32>();
        assert_relative_eq!(rslt, 0.);
        let rslt = default::<f64>();
        assert_relative_eq!(rslt, 0.);
    }

    #[test]
    fn can_autocorrelate() {
        // sin wave w/ multiple periods
        let periods = 1.;
        let points = 1000;
        let radians_per_pt = (periods * 2. * PI) / points as f64;
        let sin_wave = (0..points)
            .map(|i| (i as f64 * radians_per_pt).sin())
            .collect_vec();
        // println!("sin_wave = {:?}", sin_wave);

        let _correlations: Vec<f64> = (0..points)
            .map(|i| autocorr(sin_wave.iter(), i))
            .collect_vec();
        let correlations: Vec<f32> = (0..points)
            .map(|i| autocorr(sin_wave.iter().map(|f| *f as f32), i))
            .collect_vec();
        println!("correlations = {:?}", correlations);

        let mut fig = Figure::new();
        fig.axes2d().lines(
            correlations
                .iter()
                .enumerate()
                .map(|(i, _)| i)
                .collect_vec(),
            correlations,
            &[Caption("Sin Wave Autocorrelation")],
        );
        fig.set_pre_commands("set term dumb 100 30");
        fig.show().unwrap();
    }

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
