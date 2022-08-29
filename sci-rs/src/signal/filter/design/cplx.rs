use core::{cmp::Ordering, f64::consts::PI, ops::Mul};

use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{
    relative_degree::relative_degree, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};

/// """
/// Split into complex and real parts, combining conjugate pairs.
///
/// The 1-D input vector `z` is split up into its complex (`zc`) and real (`zr`)
/// elements. Every complex element must be part of a complex-conjugate pair,
/// which are combined into a single number (with positive imaginary part) in
/// the output. Two complex numbers are considered a conjugate pair if their
/// real and imaginary parts differ in magnitude by less than ``tol * abs(z)``.
///
/// Parameters
/// ----------
/// z : array_like
///     Vector of complex numbers to be sorted and split
/// tol : float, optional
///     Relative tolerance for testing realness and conjugate equality.
///     Default is ``100 * spacing(1)`` of `z`'s data type (i.e., 2e-14 for
///     float64)
///
/// Returns
/// -------
/// zc : ndarray
///     Complex elements of `z`, with each pair represented by a single value
///     having positive imaginary part, sorted first by real part, and then
///     by magnitude of imaginary part. The pairs are averaged when combined
///     to reduce error.
/// zr : ndarray
///     Real elements of `z` (those having imaginary part less than
///     `tol` times their magnitude), sorted by value.
///
/// Raises
/// ------
/// ValueError
///     If there are any complex numbers in `z` for which a conjugate
///     cannot be found.
///
/// See Also
/// --------
/// _cplxpair
/// """

pub fn cplxreal<F, const N: usize>(
    z: Vec<Complex<F>, N>,
    tol: Option<F>,
) -> (Vec<Complex<F>, N>, Vec<Complex<F>, N>)
where
    F: RealField + Float,
{
    if z.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Get tolerance from dtype of input
    let tol = tol.unwrap_or_else(|| F::epsilon() * F::from(100.).unwrap());

    let mut z = z;
    z.sort_unstable_by(
        |a, b| match Float::abs(a.im).partial_cmp(&Float::abs(b.im)).unwrap() {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => a.re.partial_cmp(&b.re).unwrap(),
        },
    );

    // Split reals from conjugate pairs
    let (zr, zc): (Vec<_, N>, Vec<_, N>) = z
        .iter()
        .partition(|zi| Float::abs(zi.im) <= tol * (*zi * *zi).sqrt().re);

    if zr.len() == z.len() {
        // Input is entirely real
        return (Vec::new(), zr);
    }

    // Split positive and negative halves of conjugates
    let mut zp: Vec<Complex<F>, N> = zc.iter().filter(|zi| zi.im > F::zero()).cloned().collect();
    let mut zn: Vec<Complex<F>, N> = zc.iter().filter(|zi| zi.im < F::zero()).cloned().collect();
    if zp.len() != zn.len() {
        panic!("Array contains complex value with no matchin conjugate");
    }

    // Find runs of (approximately) the same real part
    let zero_arr = [false; 1];
    let same_real: Vec<_, N> = zero_arr
        .iter()
        .cloned()
        .chain(
            zp.iter()
                .zip(zp.iter().skip(1))
                .map(|(a, b)| (a.re - b.re) <= tol * a.abs())
                .chain(zero_arr.iter().cloned()),
        )
        .collect();
    let same_real: Vec<_, N> = same_real
        .iter()
        .zip(same_real.iter().skip(1))
        .map(|(a, b)| match (a, b) {
            (true, true) => 0,
            (true, false) => -1,
            (false, true) => 1,
            (false, false) => 0,
        })
        .collect();
    let run_starts: Vec<_, N> = same_real
        .iter()
        .enumerate()
        .filter(|(i, same_real)| **same_real > 0)
        .map(|(i, _)| i)
        .collect();
    let run_stops: Vec<_, N> = same_real
        .iter()
        .enumerate()
        .filter(|(i, same_real)| **same_real < 0)
        .map(|(i, _)| i)
        .collect();

    // Sort each run by their imaginary parts
    assert!(run_starts.len() == run_stops.len());
    for i in 0..run_starts.len() {
        let start = run_starts[i];
        let stop = run_stops[i] + 1;
        let mut chunk: Vec<(Complex<F>, Complex<F>), N> = zp[start..stop]
            .iter()
            .cloned()
            .zip(zn[start..stop].iter().cloned())
            .collect();
        chunk.sort_unstable_by(|a, b| {
            match Float::abs(a.0.im).partial_cmp(&Float::abs(b.0.im)).unwrap() {
                Ordering::Less => Ordering::Less,
                Ordering::Greater => Ordering::Greater,
                Ordering::Equal => a.1.im.partial_cmp(&b.1.im).unwrap(),
            }
        });
        zp[start..stop]
            .iter_mut()
            .zip(zn[start..stop].iter_mut())
            .zip(chunk.iter())
            .for_each(|((zpi, zni), ci)| {
                *zpi = ci.0;
                *zni = ci.1;
            })
    }

    // Check that negatives match positives
    if zp
        .iter()
        .zip(zn.iter())
        .any(|(zpi, zni)| (zpi - zni.conj()).abs() > tol * zni.abs())
    {
        panic!("Array contains complex value with no matching conjugate");
    }

    // Average out numerical inaccuracy in real vs imag parts of pairs
    let zc: Vec<Complex<F>, N> = zp
        .iter()
        .zip(zn.iter())
        .map(|(zpi, zni)| (zpi + zni.conj()) / F::from(2.).unwrap())
        .collect();

    (zc, zr)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn matches_scipy_example() {}
}
