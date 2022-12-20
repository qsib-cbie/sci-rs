use core::{borrow::Borrow, iter::Sum, ops::Add};
use num_traits::Float;

// Quick select finds the `i`th smallest element with 2N comparisons
#[cfg(feature = "use_std")]
fn quickselect<B, F>(y: &Vec<B>, k: usize) -> F
where
    B: Borrow<F>,
    F: Float,
{
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
/// ```
///
#[cfg(feature = "use_std")]
pub fn median<YI, F>(y: YI) -> (F, usize)
where
    F: Float + Default,
    YI: Iterator,
    YI::Item: Borrow<F>,
{
    // Materialize the values in the iterator in order to run O(n) quick select
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
            (quickselect(&y, n / 2 - 1) + quickselect(&y, n / 2)) / F::from(2.).unwrap(),
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
/// // Flat signal perfectly correlates with itself
/// let y: [f64; 5] = [1.,2.,3.,4.,5.];
/// assert_relative_eq!(3f64, mean(y.iter()).0);
///
/// let y: &[f32] = &[];
/// assert_eq!((0f32, 0), mean(y.iter()));
///
/// ```
///
pub fn mean<YI, F>(y: YI) -> (F, usize)
where
    F: Float + Default + Add,
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

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    #[cfg(feature = "plot")]
    use gnuplot::{Figure, PlotOption::Caption};
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

        #[cfg(feature = "plot")]
        {
            let mut fig = Figure::new();
            fig.axes2d().lines(
                correlations
                    .iter()
                    .enumerate()
                    .map(|(i, _)| i)
                    .collect::<Vec<_>>(),
                correlations,
                &[Caption("Sin Wave Autocorrelation")],
            );
            fig.set_pre_commands("set term dumb 100 30");
            fig.show().unwrap();
        }
    }

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
