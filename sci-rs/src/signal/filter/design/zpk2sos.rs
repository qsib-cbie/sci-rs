use core::{cmp::Ordering, f64::consts::PI, iter::Sum, ops::Mul};

use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::{Float, Zero};

use crate::signal::filter::design::cplx::cplxreal_st;

use super::{
    zpk2tf_st, BaFormatFilter, FilterBandType, FilterOutputType, FilterType, Sos, SosFormatFilter,
    ZpkFormatFilter,
};

pub enum ZpkPairing {
    Minimal,
    Nearest,
    KeepOdd,
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
pub fn zpk2sos_st<F, const N: usize, const N2: usize>(
    zpk: ZpkFormatFilter<F, N2>,
    pairing: Option<ZpkPairing>,
    analog: Option<bool>,
) -> SosFormatFilter<F, N>
where
    F: RealField + Float + Sum,
{
    assert!(N * 2 == N2);

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
    let (zc, zr) = cplxreal_st(z, None);
    let z: Vec<Complex<F>, N2> = zc.into_iter().chain(zr.into_iter()).collect::<Vec<_, _>>();
    let (pc, pr) = cplxreal_st(p, None);
    let p: Vec<Complex<F>, N2> = pc.into_iter().chain(pr.into_iter()).collect::<Vec<_, _>>();
    let k = zpk.k;

    let idx_worst = if !analog {
        // digital: "worst" is the closest to the unit circle
        // np.argmin(np.abs(1 - np.abs(p)))
        |p: &Vec<Complex<F>, N2>| -> usize {
            p.iter()
                .enumerate()
                .map(|(i, pi)| (i, ComplexField::abs(F::one() - pi.abs())))
                .min_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(Ordering::Equal))
                .map(|(i, _)| i)
                .expect("Poles must have a min")
        }
    } else {
        //     # analog: "worst" is the closest to the imaginary axis
        //         return np.argmin(np.abs(np.real(p)))
        |p: &Vec<Complex<F>, N2>| -> usize {
            p.iter()
                .enumerate()
                .map(|(i, pi)| (i, Float::abs(pi.re)))
                .min_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(Ordering::Equal))
                .map(|(i, _)| i)
                .expect("Poles must have a min")
        }
    };

    // sos = np.zeros((n_sections, 6))
    let mut sos: Vec<Sos<F>, N> = Vec::new();

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
                let z1_idx = nearest_real_complex_idx(&z, p1, WhichNearestComplex::Real);
                let z1 = z.remove(z1_idx);
                single_zpksos(
                    Vec::from_slice(&[z1, Complex::zero()]).unwrap(),
                    Vec::from_slice(&[p1, Complex::zero()]).unwrap(),
                    F::one(),
                )
            } else if !z.is_empty() {
                let z1_idx = nearest_real_complex_idx(&z, p1, WhichNearestComplex::Real);
                let z1 = z.remove(z1_idx);
                single_zpksos(
                    Vec::from_slice(&[z1]).unwrap(),
                    Vec::from_slice(&[p1]).unwrap(),
                    F::one(),
                )
            } else {
                single_zpksos(Vec::new(), Vec::from_slice(&[p1]).unwrap(), F::one())
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
            let z1_idx = nearest_real_complex_idx(&z, p1, WhichNearestComplex::Complex);
            let z1 = z.remove(z1_idx);
            let sos_si = single_zpksos(
                Vec::from_slice(&[z1, z1.conj()]).unwrap(),
                Vec::from_slice(&[p1, p1.conj()]).unwrap(),
                F::one(),
            );
            sos.push(sos_si);
        } else {
            let p2 = if p1.im.is_zero() {
                let preal_idx = p
                    .iter()
                    .enumerate()
                    .filter(|pi| pi.1.im.is_zero())
                    .collect::<Vec<_, N2>>();
                let preal = preal_idx
                    .iter()
                    .map(|pi| pi.1)
                    .cloned()
                    .collect::<Vec<_, N2>>();
                let p2_idx = idx_worst(&preal);
                let p2_idx = preal_idx[p2_idx].0;
                drop(preal_idx);
                p.remove(p2_idx)
            } else {
                p1.conj()
            };

            // Find closest zero
            if !z.is_empty() {
                let z1_idx = nearest_real_complex_idx(&z, p1, WhichNearestComplex::Any);
                let z1 = z.remove(z1_idx);

                if !z1.im.is_zero() {
                    let sos_si = single_zpksos(
                        Vec::from_slice(&[z1, z1.conj()]).unwrap(),
                        Vec::from_slice(&[p1, p2]).unwrap(),
                        F::one(),
                    );
                    sos.push(sos_si);
                } else if !z.is_empty() {
                    let z2_idx = nearest_real_complex_idx(&z, p1, WhichNearestComplex::Real);
                    let z2 = z.remove(z2_idx);
                    assert!(z2.im.is_zero());
                    let sos_si = single_zpksos(
                        Vec::from_slice(&[z1, z2]).unwrap(),
                        Vec::from_slice(&[p1, p2]).unwrap(),
                        F::one(),
                    );
                    sos.push(sos_si);
                } else {
                    let sos_si = single_zpksos(
                        Vec::from_slice(&[z1]).unwrap(),
                        Vec::from_slice(&[p1, p2]).unwrap(),
                        F::one(),
                    );
                    sos.push(sos_si);
                }
            } else {
                let sos_si = single_zpksos(
                    Vec::from_slice(&[]).unwrap(),
                    Vec::from_slice(&[p1, p2]).unwrap(),
                    F::one(),
                );
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

// who named this stuff?
enum WhichNearestComplex {
    Real,
    Complex,
    Any,
}

/// """Get the next closest real or complex element based on distance"""
fn nearest_real_complex_idx<F, const N: usize>(
    fro: &Vec<Complex<F>, N>,
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

// """Create one second-order section from up to two zeros and poles"""
fn single_zpksos<F>(z: Vec<Complex<F>, 2>, p: Vec<Complex<F>, 2>, k: F) -> Sos<F>
where
    F: Float + RealField,
{
    let ba: BaFormatFilter<F, 3> = zpk2tf_st(&z, &p, k);
    Sos::new(ba.b, ba.a)
}
