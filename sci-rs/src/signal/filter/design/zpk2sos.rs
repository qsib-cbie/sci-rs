use core::{cmp::Ordering, f64::consts::PI, iter::Sum, ops::Mul};

use nalgebra::{Complex, ComplexField, RealField};
use num_traits::{Float, Zero};

#[cfg(feature = "alloc")]
use crate::signal::filter::design::cplx::cplxreal_dyn;

#[cfg(feature = "alloc")]
use super::{
    zpk2tf_dyn, BaFormatFilter, FilterBandType, FilterOutputType, FilterType, Sos, SosFormatFilter,
    ZpkFormatFilter,
};

#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

///
/// Choice of zeros/poles pairing for the zpk2sos conversion
///
/// Matches scipy.signal.zpk2sos pairing
///
pub enum ZpkPairing {
    /// Default
    Minimal,
    /// Nearest
    Nearest,
}

enum WhichNearestComplex {
    Real,
    Complex,
    Any,
}

/// Return second-order sections from zeros, poles, and gain of a system
///
/// Returns
/// -------
/// sos : ndarray
///     Array of second-order filter coefficients, with shape
///     ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
///     specification.
///
///

#[cfg(feature = "alloc")]
pub fn zpk2sos_dyn<F>(
    order: usize,
    zpk: ZpkFormatFilter<F>,
    pairing: Option<ZpkPairing>,
    analog: Option<bool>,
) -> SosFormatFilter<F>
where
    F: RealField + Float + Sum,
{
    let analog = analog.unwrap_or(false);
    let pairing = pairing.unwrap_or(if analog {
        ZpkPairing::Minimal
    } else {
        ZpkPairing::Nearest
    });

    if analog && !matches!(pairing, ZpkPairing::Minimal) {
        panic!("for analog zpk2sos conversion, pairing must be minimal");
    }

    if zpk.z.is_empty() && zpk.p.is_empty() {
        if !analog {
            return SosFormatFilter {
                sos: Vec::from_iter(
                    [Sos::new(
                        [zpk.k, F::zero(), F::zero()],
                        [F::one(), F::zero(), F::zero()],
                    )]
                    .iter()
                    .cloned(),
                ),
            };
        } else {
            return SosFormatFilter {
                sos: Vec::from_iter(
                    [Sos::new(
                        [F::zero(), F::zero(), zpk.k],
                        [F::zero(), F::zero(), F::one()],
                    )]
                    .iter()
                    .cloned(),
                ),
            };
        }
    }

    let mut z = zpk.z;
    let mut p = zpk.p;
    let n_sections = if !matches!(pairing, ZpkPairing::Minimal) {
        // ensure we have the same number of poles and zeros, and make copies
        let z_p = z.len() as isize - p.len() as isize;
        if z_p > 0 {
            p.extend((0..z_p).map(|_| Complex::zero()));
        }
        let p_z = p.len() as isize - z.len() as isize;
        if p_z > 0 {
            z.extend((0..p_z).map(|_| Complex::zero()));
        }

        // TODO: Why is this logic funky?
        let n_sections = (p.len().max(z.len()) + 1) / 2;

        if p.len() % 2 == 1 && matches!(pairing, ZpkPairing::Nearest) {
            z.push(Complex::zero());
            p.push(Complex::zero());
        }

        assert!(z.len() == p.len());

        n_sections
    } else {
        if p.len() < z.len() {
            panic!("for analog zpk2sos conversion, must have len(p)>=len(z)");
        }
        (p.len() + 1) / 2
    };

    // Ensure we have complex conjugate pairs
    // (note that _cplxreal only gives us one element of each complex pair):
    let (zc, zr) = cplxreal_dyn(z, None);
    let z: Vec<Complex<F>> = zc.into_iter().chain(zr).collect::<Vec<_>>();
    let (pc, pr) = cplxreal_dyn(p, None);
    let p: Vec<Complex<F>> = pc.into_iter().chain(pr).collect::<Vec<_>>();
    let k = zpk.k;

    let idx_worst = if !analog {
        // digital: "worst" is the closest to the unit circle
        // np.argmin(np.abs(1 - np.abs(p)))
        |p: &Vec<Complex<F>>| -> usize {
            p.iter()
                .enumerate()
                .map(|(i, pi)| (i, ComplexField::abs(F::one() - pi.abs())))
                .min_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(Ordering::Equal))
                .map(|(i, _)| i)
                .expect("Poles must have a min")
        }
    } else {
        // analog: "worst" is the closest to the imaginary axis
        // np.argmin(np.abs(np.real(p)))
        |p: &Vec<Complex<F>>| -> usize {
            p.iter()
                .enumerate()
                .map(|(i, pi)| (i, Float::abs(pi.re)))
                .min_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(Ordering::Equal))
                .map(|(i, _)| i)
                .expect("Poles must have a min")
        }
    };

    // sos = np.zeros((n_sections, 6))
    let mut sos: Vec<Sos<F>> = Vec::new();

    // Construct the system s..t the "worst" are first
    let mut z = z;
    let mut p = p;
    for si in 0..n_sections {
        // Select the next "worst" pole
        let p1_idx = idx_worst(&p);
        let p1 = p.remove(p1_idx);

        // Pair that pole with a zero
        if p1.im.is_zero() && p.iter().map(|pi| pi.re).sum::<F>().is_zero() {
            // Special case (1): last remaining real pole
            let sos_si = if matches!(pairing, ZpkPairing::Minimal) {
                let z1_idx = nearest_real_complex_idx_dyn(&z, p1, WhichNearestComplex::Real);
                let z1 = z.remove(z1_idx);
                single_zpksos_dyn(
                    vec![z1, Complex::zero()],
                    vec![p1, Complex::zero()],
                    F::one(),
                )
            } else if !z.is_empty() {
                let z1_idx = nearest_real_complex_idx_dyn(&z, p1, WhichNearestComplex::Real);
                let z1 = z.remove(z1_idx);
                single_zpksos_dyn(vec![z1], vec![p1], F::one())
            } else {
                single_zpksos_dyn(Vec::new(), vec![p1], F::one())
            };
            sos.push(sos_si);
        } else if p.len() + 1 == z.len()
            && !p1.im.is_zero()
            && p.iter().filter(|pi| pi.im.is_zero()).count() == 1
            && z.iter().filter(|zi| zi.im.is_zero()).count() == 1
        {
            // Special case (2): there's one real pole and one real zero
            // left, and an equal number of poles and zeros to pair up.
            // We *must* pair with a complex zero
            let z1_idx = nearest_real_complex_idx_dyn(&z, p1, WhichNearestComplex::Complex);
            let z1 = z.remove(z1_idx);
            let sos_si = single_zpksos_dyn(vec![z1, z1.conj()], vec![p1, p1.conj()], F::one());
            sos.push(sos_si);
        } else {
            let p2 = if p1.im.is_zero() {
                let preal_idx = p
                    .iter()
                    .enumerate()
                    .filter(|pi| pi.1.im.is_zero())
                    .collect::<Vec<_>>();
                let preal = preal_idx.iter().map(|pi| pi.1).cloned().collect::<Vec<_>>();
                let p2_idx = idx_worst(&preal);
                let p2_idx = preal_idx[p2_idx].0;
                drop(preal_idx);
                p.remove(p2_idx)
            } else {
                p1.conj()
            };

            // Find closest zero
            if !z.is_empty() {
                let z1_idx = nearest_real_complex_idx_dyn(&z, p1, WhichNearestComplex::Any);
                let z1 = z.remove(z1_idx);

                if !z1.im.is_zero() {
                    let sos_si = single_zpksos_dyn(vec![z1, z1.conj()], vec![p1, p2], F::one());
                    sos.push(sos_si);
                } else if !z.is_empty() {
                    let z2_idx = nearest_real_complex_idx_dyn(&z, p1, WhichNearestComplex::Real);
                    let z2 = z.remove(z2_idx);
                    assert!(z2.im.is_zero());
                    let sos_si = single_zpksos_dyn(vec![z1, z2], vec![p1, p2], F::one());
                    sos.push(sos_si);
                } else {
                    let sos_si = single_zpksos_dyn(vec![z1], vec![p1, p2], F::one());
                    sos.push(sos_si);
                }
            } else {
                let sos_si = single_zpksos_dyn(vec![], vec![p1, p2], F::one());
                sos.push(sos_si);
            }
        }
    }

    // Reverse so the "worst" are last
    sos.reverse();

    assert!(p.len() == z.len());
    assert!(p.is_empty());

    // Put the gain in the first sos
    for bi in sos[0].b.iter_mut() {
        *bi *= zpk.k;
    }

    SosFormatFilter { sos }
}

/// """Get the next closest real or complex element based on distance"""
#[cfg(feature = "alloc")]
fn nearest_real_complex_idx_dyn<F>(
    fro: &[Complex<F>],
    to: Complex<F>,
    which: WhichNearestComplex,
) -> usize
where
    F: Float + RealField,
{
    let order = fro.iter().map(|fi| (Complex::abs(fi - to), fi)).enumerate();
    match which {
        WhichNearestComplex::Real => order
            .filter(|ai| ai.1 .1.im.is_zero())
            .min_by(|a, b| a.1 .0.partial_cmp(&b.1 .0).unwrap_or(Ordering::Equal))
            .map(|a| a.0)
            .expect("Min must exist"),
        WhichNearestComplex::Complex => order
            .filter(|ai| !ai.1 .1.im.is_zero())
            .min_by(|a, b| a.1 .0.partial_cmp(&b.1 .0).unwrap_or(Ordering::Equal))
            .map(|a| a.0)
            .expect("Min must exist"),
        WhichNearestComplex::Any => order
            .min_by(|a, b| a.1 .0.partial_cmp(&b.1 .0).unwrap_or(Ordering::Equal))
            .map(|a| a.0)
            .expect("Min must exist"),
    }
}

/// """Create one second-order section from up to two zeros and poles"""
#[cfg(feature = "alloc")]
fn single_zpksos_dyn<F>(z: Vec<Complex<F>>, p: Vec<Complex<F>>, k: F) -> Sos<F>
where
    F: Float + RealField,
{
    let ba: BaFormatFilter<F> = zpk2tf_dyn(2, &z, &p, k);
    if ba.b.len() != 3 || ba.a.len() != 3 {
        panic!(
            "SOS must have 3 coefficients has {} and {}",
            ba.b.len(),
            ba.a.len()
        );
    }
    Sos::new([ba.b[0], ba.b[1], ba.b[2]], [ba.a[0], ba.a[1], ba.a[2]])
}
