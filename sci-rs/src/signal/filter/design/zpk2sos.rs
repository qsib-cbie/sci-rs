use core::{f64::consts::PI, ops::Mul};

use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{FilterBandType, FilterOutputType, FilterType, Sos, SosFormatFilter, ZpkFormatFilter};

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
    analog: bool,
) -> SosFormatFilter<F, N>
where
    F: RealField + Float,
{
    assert!(N * 2 == N2);
    todo!()
    // # TODO in the near future:
    // # 1. Add SOS capability to `filtfilt`, `freqz`, etc. somehow (#3259).
    // # 2. Make `decimate` use `sosfilt` instead of `lfilter`.
    // # 3. Make sosfilt automatically simplify sections to first order
    // #    when possible. Note this might make `sosfiltfilt` a bit harder (ICs).
    // # 4. Further optimizations of the section ordering / pole-zero pairing.
    // # See the wiki for other potential issues.

    // if pairing is None:
    //     pairing = 'minimal' if analog else 'nearest'

    // valid_pairings = ['nearest', 'keep_odd', 'minimal']
    // if pairing not in valid_pairings:
    //     raise ValueError('pairing must be one of %s, not %s'
    //                      % (valid_pairings, pairing))

    // if analog and pairing != 'minimal':
    //     raise ValueError('for analog zpk2sos conversion, '
    //                      'pairing must be "minimal"')

    // if len(z) == len(p) == 0:
    //     if not analog:
    //         return np.array([[k, 0., 0., 1., 0., 0.]])
    //     else:
    //         return np.array([[0., 0., k, 0., 0., 1.]])

    // if pairing != 'minimal':
    //     # ensure we have the same number of poles and zeros, and make copies
    //     p = np.concatenate((p, np.zeros(max(len(z) - len(p), 0))))
    //     z = np.concatenate((z, np.zeros(max(len(p) - len(z), 0))))
    //     n_sections = (max(len(p), len(z)) + 1) // 2

    //     if len(p) % 2 == 1 and pairing == 'nearest':
    //         p = np.concatenate((p, [0.]))
    //         z = np.concatenate((z, [0.]))
    //     assert len(p) == len(z)
    // else:
    //     if len(p) < len(z):
    //         raise ValueError('for analog zpk2sos conversion, '
    //                          'must have len(p)>=len(z)')

    //     n_sections = (len(p) + 1) // 2

    // # Ensure we have complex conjugate pairs
    // # (note that _cplxreal only gives us one element of each complex pair):
    // z = np.concatenate(_cplxreal(z))
    // p = np.concatenate(_cplxreal(p))
    // if not np.isreal(k):
    //     raise ValueError('k must be real')
    // k = k.real

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
}
