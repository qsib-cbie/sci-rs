use core::{f64::consts::PI, ops::Mul};

use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::{Float, Zero};

use crate::signal::filter::design::cplx::cplxreal;

use super::{FilterBandType, FilterOutputType, FilterType, Sos, SosFormatFilter, ZpkFormatFilter};

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
pub fn zpk2sos<F, const N: usize, const N2: usize>(
    zpk: ZpkFormatFilter<F, N2>,
    pairing: Option<ZpkPairing>,
    analog: Option<bool>,
) -> SosFormatFilter<F, N>
where
    F: RealField + Float,
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

    if zpk.z.len() == 0 && zpk.p.len() == 0 {
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
    let (zc, zr) = cplxreal(z, None);
    let z: Vec<Complex<F>, N> = zc.into_iter().chain(zr.into_iter()).collect::<Vec<_, _>>();
    let (pc, pr) = cplxreal(p, None);
    let p: Vec<Complex<F>, N> = pc.into_iter().chain(pr.into_iter()).collect::<Vec<_, _>>();
    let k = zpk.k;

    // if not analog:
    //     # digital: "worst" is the closest to the unit circle
    //     def idx_worst(p):
    //         return np.argmin(np.abs(1 - np.abs(p)))
    // else:
    //     # analog: "worst" is the closest to the imaginary axis
    //     def idx_worst(p):
    //         return np.argmin(np.abs(np.real(p)))

    // sos = np.zeros((n_sections, 6))

    // # Construct the system, reversing order so the "worst" are last
    // for si in range(n_sections-1, -1, -1):
    //     # Select the next "worst" pole
    //     p1_idx = idx_worst(p)
    //     p1 = p[p1_idx]
    //     p = np.delete(p, p1_idx)

    //     # Pair that pole with a zero

    //     if np.isreal(p1) and np.isreal(p).sum() == 0:
    //         # Special case (1): last remaining real pole
    //         if pairing != 'minimal':
    //             z1_idx = _nearest_real_complex_idx(z, p1, 'real')
    //             z1 = z[z1_idx]
    //             z = np.delete(z, z1_idx)
    //             sos[si] = _single_zpksos([z1, 0], [p1, 0], 1)
    //         elif len(z) > 0:
    //             z1_idx = _nearest_real_complex_idx(z, p1, 'real')
    //             z1 = z[z1_idx]
    //             z = np.delete(z, z1_idx)
    //             sos[si] = _single_zpksos([z1], [p1], 1)
    //         else:
    //             sos[si] = _single_zpksos([], [p1], 1)

    //     elif (len(p) + 1 == len(z)
    //           and not np.isreal(p1)
    //           and np.isreal(p).sum() == 1
    //           and np.isreal(z).sum() == 1):

    //         # Special case (2): there's one real pole and one real zero
    //         # left, and an equal number of poles and zeros to pair up.
    //         # We *must* pair with a complex zero

    //         z1_idx = _nearest_real_complex_idx(z, p1, 'complex')
    //         z1 = z[z1_idx]
    //         z = np.delete(z, z1_idx)
    //         sos[si] = _single_zpksos([z1, z1.conj()], [p1, p1.conj()], 1)

    //     else:
    //         if np.isreal(p1):
    //             prealidx = np.flatnonzero(np.isreal(p))
    //             p2_idx = prealidx[idx_worst(p[prealidx])]
    //             p2 = p[p2_idx]
    //             p = np.delete(p, p2_idx)
    //         else:
    //             p2 = p1.conj()

    //         # find closest zero
    //         if len(z) > 0:
    //             z1_idx = _nearest_real_complex_idx(z, p1, 'any')
    //             z1 = z[z1_idx]
    //             z = np.delete(z, z1_idx)

    //             if not np.isreal(z1):
    //                 sos[si] = _single_zpksos([z1, z1.conj()], [p1, p2], 1)
    //             else:
    //                 if len(z) > 0:
    //                     z2_idx = _nearest_real_complex_idx(z, p1, 'real')
    //                     z2 = z[z2_idx]
    //                     assert np.isreal(z2)
    //                     z = np.delete(z, z2_idx)
    //                     sos[si] = _single_zpksos([z1, z2], [p1, p2], 1)
    //                 else:
    //                     sos[si] = _single_zpksos([z1], [p1, p2], 1)
    //         else:
    //             # no more zeros
    //             sos[si] = _single_zpksos([], [p1, p2], 1)

    // assert len(p) == len(z) == 0  # we've consumed all poles and zeros
    // del p, z

    // # put gain in first sos
    // sos[0][:3] *= k
    // return sos
    todo!()
}
