use core::{borrow::Borrow, iter::Sum, ops::Add};
use itertools::Itertools;
use num_traits::{Float, Num, NumCast};

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

// Quick select finds the `i`th smallest element with 2N comparisons
#[cfg(feature = "alloc")]
fn quickselect<B, T>(y: &Vec<B>, k: usize) -> T
where
    B: Borrow<T>,
    T: Num + NumCast + PartialOrd + Copy,
{
    use num_traits::{Num, NumCast};

    let n = y.len();
    if n == 1 {
        return *y[0].borrow();
    }

    let pivot = y.get(n / 2).unwrap().borrow();
    let lower = y
        .iter()
        .filter(|yi| *(*yi).borrow() < *pivot)
        .map(|yi| *yi.borrow())
        .collect::<Vec<_>>();
    let lowers = lower.len();
    let upper = y
        .iter()
        .filter(|yi| *(*yi).borrow() > *pivot)
        .map(|yi| *yi.borrow())
        .collect::<Vec<_>>();
    let uppers = upper.len();
    let pivots = n - lowers - uppers;

    if k < lowers {
        quickselect(&lower, k)
    } else if k < lowers + pivots {
        *pivot
    } else {
        quickselect(&upper, k - lowers - pivots)
    }
}

///
/// Compute the median of the signal, `y`
///
/// Return the median and the number of points averaged
///
/// ```
/// use approx::assert_relative_eq;
/// use sci_rs::stats::median;
///
/// let y: [f64; 5] = [1.,2.,3.,4.,5.];
/// assert_relative_eq!(3f64, median(y.iter()).0);
///
/// let y: [i32; 5] = [1,2,3,4,5];
/// assert_eq!(3, median(y.iter()).0);
///
/// let y: [f64; 4] = [1.,2.,3.,4.];
/// assert_relative_eq!(2.5f64, median(y.iter()).0);
///
/// let y: [f32; 5] = [3.,1.,4.,2.,5.];
/// assert_relative_eq!(3f32, median(y.iter()).0);
///
/// let y: [f64; 6] = [3.,1.,4.,2.,3.,5.];
/// assert_relative_eq!(3f64, median(y.iter()).0);
///
/// let y: &[f32] = &[];
/// assert_eq!((0f32, 0), median(y.iter()));
///
/// let y: &[f32] = &[1.];
/// assert_eq!((1f32, 1), median(y.iter()));
///
/// let y: [i64; 4] = [1,2,3,4];
/// assert_eq!(2i64, median(y.iter()).0);
///
/// ```
///
#[cfg(feature = "alloc")]
pub fn median<YI, T>(y: YI) -> (T, usize)
where
    T: Num + NumCast + PartialOrd + Copy + Default,
    YI: Iterator,
    YI::Item: Borrow<T>,
{
    // Materialize the values in the iterator in order to run O(n) quick select

    use num_traits::NumCast;
    let y = y.collect::<Vec<_>>();
    let n = y.len();

    if n == 0 {
        Default::default()
    } else if n == 1 {
        (*y[0].borrow(), 1)
    } else if n % 2 == 1 {
        (quickselect(&y, n / 2), n)
    } else {
        (
            (quickselect(&y, n / 2 - 1) + quickselect(&y, n / 2)) / T::from(2).unwrap(),
            n,
        )
    }
}

///
/// Compute the mean of the signal, `y`
///
/// Return the mean and the number of points averaged
///
/// ```
/// use approx::assert_relative_eq;
/// use sci_rs::stats::mean;
///
/// let y: [f64; 5] = [1.,2.,3.,4.,5.];
/// assert_relative_eq!(3f64, mean(y.iter()).0);
///
/// let y: [i64; 5] = [1,2,3,4,5];
/// assert_eq!(3i64, mean(y.iter()).0);
///
/// let y: &[f32] = &[];
/// assert_eq!((0f32, 0), mean(y.iter()));
///
/// ```
///
pub fn mean<YI, F>(y: YI) -> (F, usize)
where
    F: Num + NumCast + Default + Copy + Add,
    YI: Iterator,
    YI::Item: Borrow<F>,
{
    let (sum, count) = y.fold(Default::default(), |acc: (F, usize), yi| {
        (acc.0 + *yi.borrow(), acc.1 + 1)
    });
    if count > 0 {
        (sum / F::from(count).unwrap(), count)
    } else {
        Default::default()
    }
}

///
/// Compute the variance of the signal, `y`
///
/// Return the variance and the number of points averaged
///
/// ```
/// use approx::assert_relative_eq;
/// use sci_rs::stats::variance;
///
/// let y: [f64; 5] = [1.,2.,3.,4.,5.];
/// assert_relative_eq!(2f64, variance(y.iter()).0);
///
/// let y: &[f32] = &[];
/// assert_eq!((0f32, 0), variance(y.iter()));
///
/// ```
///
pub fn variance<YI, F>(y: YI) -> (F, usize)
where
    F: Float + Default + Sum,
    YI: Iterator + Clone,
    YI::Item: Borrow<F>,
{
    let (avg, n) = mean(y.clone());
    let sum: F = y
        .map(|f| {
            let delta = *f.borrow() - avg;
            delta * delta
        })
        .sum::<F>();
    if n > 0 {
        (sum / F::from(n).unwrap(), n)
    } else {
        Default::default()
    }
}

///
/// Compute the standard deviation of the signal, `y`
///
/// Return the standard deviation and the number of points averaged
///
/// ```
/// use approx::assert_relative_eq;
/// use sci_rs::stats::stdev;
///
/// let y: [f64; 5] = [1.,2.,3.,4.,5.];
/// assert_relative_eq!(1.41421356237, stdev(y.iter()).0, max_relative = 1e-8);
///
/// let y: &[f32] = &[];
/// assert_eq!((0f32, 0), stdev(y.iter()));
///
/// ```
pub fn stdev<YI, F>(y: YI) -> (F, usize)
where
    F: Float + Default + Sum,
    YI: Iterator + Clone,
    YI::Item: Borrow<F>,
{
    match variance(y) {
        (_, 0) => Default::default(),
        (v, n) => (v.sqrt(), n),
    }
}

///
/// Autocorrelate the signal `y` with lag `k`,
/// using 1/N formulation
///
/// <https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm>
///
pub fn autocorr<YI, F>(y: YI, k: usize) -> F
where
    F: Float + Add + Sum + Default,
    YI: Iterator + Clone,
    YI::Item: Borrow<F>,
{
    let (avg, n) = mean(y.clone());
    let n = F::from(n).unwrap();
    let (var, _) = variance(y.clone());
    let autocovariance: F = y
        .clone()
        .zip(y.skip(k))
        .map(|(fi, fik)| (*fi.borrow() - avg) * (*fik.borrow() - avg))
        .sum::<F>()
        / n;
    autocovariance / var
}

///
/// Unscaled tiled autocorrelation of signal `y` with itself into `x`.
///
/// This skips variance normalization and only computes lags in `SKIP..SKIP+x.len()`
///
/// The autocorrelation is not normalized by 1/y.len() or variance. The variance of the signal
/// is returned. The returned variance may be used to normalize lags of interest after the fact.
///
/// <https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm>
///
pub fn autocorr_fast32<const N: usize, const M: usize, const SKIP: usize>(
    y: &mut [f32; N],
    x: &mut [f32; M],
) -> f32 {
    assert!(N >= M + SKIP);

    // Subtract the mean
    let sum = y.iter().sum::<f32>();
    let avg = sum / y.len() as f32;
    y.iter_mut().for_each(|yi| *yi = (*yi - avg));

    // Compute the variance of the signal
    let var = y.iter().map(|yi| yi * yi).sum::<f32>() / y.len() as f32;

    // Compute the autocorrelation for lag 1 to lag n
    let lag_skip = y.len() - x.len();
    for (h, xi) in (SKIP..y.len()).zip(x.iter_mut()) {
        let left = &y[..y.len() - h];
        let right = &y[h..];
        const TILE: usize = 4;
        let left = left.chunks_exact(TILE);
        let right = right.chunks_exact(TILE);
        *xi = left
            .remainder()
            .iter()
            .zip(right.remainder().iter())
            .map(|(a, b)| a * b)
            .sum::<f32>();
        *xi = left
            .zip(right)
            .map(|(left, right)| {
                left.iter()
                    .zip(right.iter())
                    .map(|(a, b)| a * b)
                    .sum::<f32>()
            })
            .sum();
    }

    var
}

///
/// Root Mean Square (RMS) of signal `y`.
///
/// It is assumed that the mean of the signal is zero.
///
pub fn rms_fast32<const N: usize>(y: &[f32; N]) -> f32 {
    const TILE: usize = 4;
    let tiles = y.chunks_exact(TILE);
    let sum = tiles.remainder().iter().map(|yi| yi * yi).sum::<f32>()
        + tiles
            .map(|yi| yi.iter().map(|yi| yi * yi).sum::<f32>())
            .sum::<f32>();
    (sum / y.len() as f32).sqrt()
}

///
/// Produce an iterator yielding the lag difference, yi1 - yi0,
///
/// <https://www.itl.nist.gov/div898/handbook/eda/section3/lagplot.htm>
///
/// ```
/// use approx::assert_relative_eq;
/// use sci_rs::stats::lag_diff;
///
/// // Flat signal perfectly correlates with itself
/// let y: [f64; 4] = [1.,2.,4.,7.];
/// let z = lag_diff(y.iter()).collect::<Vec<_>>();
/// for i in 0..3 {
///     assert_relative_eq!(i as f64 + 1f64, z[i]);
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

///
/// Compute the root mean square of successive differences
///
/// ```
/// use approx::assert_relative_eq;
/// use sci_rs::stats::rmssd;
///
/// // Differences are 1, 2, 3
/// // Square differences are 1, 4, 9
/// // Mean is 4.666666666666667
/// // RMSSD is 2.1602468995
/// let y: [f64; 4] = [1.,2.,4.,7.];
/// assert_relative_eq!(2.1602468995, rmssd(y.iter()), max_relative = 1e-8);
/// ```
///
pub fn rmssd<YI, F>(y: YI) -> F
where
    F: Float + Add + Sum + Default,
    YI: Iterator + Clone,
    YI::Item: Borrow<F> + Copy,
{
    let square_diffs = y
        .tuple_windows()
        .map(|(yi0, yi1)| (*yi1.borrow() - *yi0.borrow()).powi(2));
    let (sum, n): (F, usize) = mean(square_diffs);
    sum.sqrt()
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use std::f64::consts::PI;
    use std::vec::Vec;

    use super::*;

    #[test]
    fn can_median() {
        let y: [f64; 4] = [1., 2., 3., 4.];
        println!("y = {:?}", y);
        println!("y = {:?}", median::<_, f64>(y.iter()));
        assert_relative_eq!(2.5, median::<_, f64>(y.iter()).0);
        let y: [f64; 5] = [1., 2., 3., 4., 5.];
        println!("y = {:?}", y);
        println!("y = {:?}", median::<_, f64>(y.iter()));
        assert_relative_eq!(3.0, median::<_, f64>(y.iter()).0);
    }

    #[test]
    fn can_autocorrelate() {
        // sin wave w/ multiple periods
        let periods = 1.;
        let points = 100;
        let radians_per_pt = (periods * 2. * PI) / points as f64;
        let sin_wave = (0..points)
            .map(|i| (i as f64 * radians_per_pt).sin())
            .collect::<Vec<_>>();
        // println!("sin_wave = {:?}", sin_wave);

        let _correlations: Vec<f64> = (0..points)
            .map(|i| autocorr(sin_wave.iter(), i))
            .collect::<Vec<_>>();
        let correlations: Vec<f32> = (0..points)
            .map(|i| autocorr(sin_wave.iter().map(|f| *f as f32), i))
            .collect::<Vec<_>>();
        println!("correlations = {:?}", correlations);
    }

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
