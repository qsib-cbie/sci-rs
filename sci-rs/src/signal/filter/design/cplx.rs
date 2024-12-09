use core::{cmp::Ordering, f64::consts::PI, ops::Mul};
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

#[cfg(feature = "alloc")]
use super::{
    relative_degree::relative_degree_dyn, FilterBandType, FilterOutputType, FilterType, Sos,
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
#[cfg(feature = "alloc")]
pub fn cplxreal_dyn<F>(z: Vec<Complex<F>>, tol: Option<F>) -> (Vec<Complex<F>>, Vec<Complex<F>>)
where
    F: RealField + Float,
{
    if z.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Get tolerance from dtype of input
    let tol = tol.unwrap_or_else(|| F::epsilon() * F::from(100.).unwrap());

    let mut z = z;
    z.sort_unstable_by(|a, b| match a.re.partial_cmp(&b.re).unwrap() {
        Ordering::Less => Ordering::Less,
        Ordering::Greater => Ordering::Greater,
        Ordering::Equal => Float::abs(a.im).partial_cmp(&Float::abs(b.im)).unwrap(),
    });

    // Split reals from conjugate pairs
    let (zr, zc): (Vec<_>, Vec<_>) = z
        .iter()
        .partition(|zi| Float::abs(zi.im) <= tol * (*zi * *zi).sqrt().re);

    if zr.len() == z.len() {
        // Input is entirely real
        return (Vec::new(), zr);
    }

    // Split positive and negative halves of conjugates
    let mut zp: Vec<Complex<F>> = zc.iter().filter(|zi| zi.im > F::zero()).cloned().collect();
    let mut zn: Vec<Complex<F>> = zc.iter().filter(|zi| zi.im < F::zero()).cloned().collect();
    if zp.len() != zn.len() {
        panic!("Array contains complex value with no matching conjugate");
    }

    // Find runs of (approximately) the same real part
    let zero_arr = [false; 1];
    let same_real: Vec<_> = zero_arr
        .iter()
        .cloned()
        .chain(
            zp.iter()
                .zip(zp.iter().skip(1))
                .map(|(a, b)| (b.re - a.re) <= tol * a.abs())
                .chain(zero_arr.iter().cloned()),
        )
        .collect();
    let same_real: Vec<_> = same_real
        .iter()
        .zip(same_real.iter().skip(1))
        .map(|(a, b)| match (a, b) {
            (true, true) => 0,
            (true, false) => -1,
            (false, true) => 1,
            (false, false) => 0,
        })
        .collect();
    let run_starts: Vec<_> = same_real
        .iter()
        .enumerate()
        .filter(|(i, same_real)| **same_real > 0)
        .map(|(i, _)| i)
        .collect();
    let run_stops: Vec<_> = same_real
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
        let mut chunk: Vec<(Complex<F>, Complex<F>)> = zp[start..stop]
            .iter()
            .cloned()
            .zip(zn[start..stop].iter().cloned())
            .collect();
        chunk.sort_unstable_by(|a, b| {
            match Float::abs(a.1.im).partial_cmp(&Float::abs(b.1.im)).unwrap() {
                Ordering::Less => Ordering::Less,
                Ordering::Greater => Ordering::Greater,
                Ordering::Equal => a.0.im.partial_cmp(&b.0.im).unwrap(),
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
    let zc: Vec<Complex<F>> = zp
        .iter()
        .zip(zn.iter())
        .map(|(zpi, zni)| (zpi + zni.conj()) / F::from(2.).unwrap())
        .collect();

    (zc, zr)
}

#[cfg(feature = "alloc")]
pub fn sort_cplx_dyn<F: ComplexField>(x: &mut [F]) {
    x.sort_unstable_by(|a, b| {
        match a
            .clone()
            .real()
            .partial_cmp(&b.clone().real())
            .expect("Reals must be orderable")
        {
            Ordering::Equal => a
                .clone()
                .imaginary()
                .partial_cmp(&b.clone().imaginary())
                .expect("Imaginaries must be orderable"),
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
        }
    });
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[cfg(all(feature = "alloc", feature = "std"))]
    #[test]
    fn matches_scipy_example() {
        let z = vec![
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(-1., 0.),
            Complex::new(-1., 0.),
            Complex::new(-1., 0.),
            Complex::new(-1., 0.),
        ];

        let expected_zc: Vec<Complex<f64>> = Vec::new();
        let expected_zr = vec![
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
        ];

        let (zc, zr) = cplxreal_dyn(z, None);

        assert_eq!(expected_zc.len(), zc.len());
        expected_zc.iter().zip(zc.iter()).for_each(|(e, a)| {
            assert_relative_eq!(e.re, a.re, max_relative = 1e-7);
            assert_relative_eq!(e.im, a.im, max_relative = 1e-7);
        });

        let z = vec![
            Complex::new(0.98924866, -0.03710237),
            Complex::new(0.96189799, -0.03364097),
            Complex::new(0.96189799, 0.03364097),
            Complex::new(0.98924866, 0.03710237),
            Complex::new(0.93873849, 0.16792939),
            Complex::new(0.89956011, 0.08396115),
            Complex::new(0.89956011, -0.08396115),
            Complex::new(0.93873849, -0.16792939),
        ];

        let expected_zc = vec![
            Complex::new(0.89956011, 0.08396115),
            Complex::new(0.93873849, 0.16792939),
            Complex::new(0.96189799, 0.03364097),
            Complex::new(0.98924866, 0.03710237),
        ];
        let expected_zr: Vec<Complex<f64>> = Vec::new();

        let (zc, zr) = cplxreal_dyn(z, None);

        assert_eq!(expected_zc.len(), zc.len());
        expected_zc.iter().zip(zc.iter()).for_each(|(e, a)| {
            assert_relative_eq!(e.re, a.re, max_relative = 1e-7);
            assert_relative_eq!(e.im, a.im, max_relative = 1e-7);
        });

        // TODO: add test cases with > 0 run_starts and run_stops
    }
}
