use core::{f64::consts::PI, ops::Mul};

use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{FilterBandType, FilterOutput, FilterOutputType, FilterType, Sos, Zpk};

/// """Return second-order sections from zeros, poles, and gain of a system
///
/// Parameters
/// ----------
/// z : array_like
///     Zeros of the transfer function.
/// p : array_like
///     Poles of the transfer function.
/// k : float
///     System gain.
/// pairing : {None, 'nearest', 'keep_odd', 'minimal'}, optional
///     The method to use to combine pairs of poles and zeros into sections.
///     If analog is False and pairing is None, pairing is set to 'nearest';
///     if analog is True, pairing must be 'minimal', and is set to that if
///     it is None.
/// analog : bool, optional
///     If True, system is analog, otherwise discrete.
///
///     .. versionadded:: 1.8.0
///
/// Returns
/// -------
/// sos : ndarray
///     Array of second-order filter coefficients, with shape
///     ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
///     specification.
///
/// See Also
/// --------
/// sosfilt
///
/// Notes
/// -----
/// The algorithm used to convert ZPK to SOS format is designed to
/// minimize errors due to numerical precision issues. The pairing
/// algorithm attempts to minimize the peak gain of each biquadratic
/// section. This is done by pairing poles with the nearest zeros, starting
/// with the poles closest to the unit circle for discrete-time systems, and
/// poles closest to the imaginary axis for continuous-time systems.
///
/// ``pairing='minimal'`` outputs may not be suitable for `sosfilt`,
/// and ``analog=True`` outputs will never be suitable for `sosfilt`.
///
/// *Algorithms*
///
/// The steps in the ``pairing='nearest'``, ``pairing='keep_odd'``,
/// and ``pairing='minimal'`` algorithms are mostly shared. The
/// ``'nearest'`` algorithm attempts to minimize the peak gain, while
/// ``'keep_odd'`` minimizes peak gain under the constraint that
/// odd-order systems should retain one section as first order.
/// ``'minimal'`` is similar to ``'keep_odd'``, but no additional
/// poles or zeros are introduced
///
/// The algorithm steps are as follows:
///
/// As a pre-processing step for ``pairing='nearest'``,
/// ``pairing='keep_odd'``, add poles or zeros to the origin as
/// necessary to obtain the same number of poles and zeros for
/// pairing.  If ``pairing == 'nearest'`` and there are an odd number
/// of poles, add an additional pole and a zero at the origin.
///
/// The following steps are then iterated over until no more poles or
/// zeros remain:
///
/// 1. Take the (next remaining) pole (complex or real) closest to the
///    unit circle (or imaginary axis, for ``analog=True``) to
///    begin a new filter section.
///
/// 2. If the pole is real and there are no other remaining real poles [#]_,
///    add the closest real zero to the section and leave it as a first
///    order section. Note that after this step we are guaranteed to be
///    left with an even number of real poles, complex poles, real zeros,
///    and complex zeros for subsequent pairing iterations.
///
/// 3. Else:
///
///     1. If the pole is complex and the zero is the only remaining real
///        zero*, then pair the pole with the *next* closest zero
///        (guaranteed to be complex). This is necessary to ensure that
///        there will be a real zero remaining to eventually create a
///        first-order section (thus keeping the odd order).
///
///     2. Else pair the pole with the closest remaining zero (complex or
///        real).
///
///     3. Proceed to complete the second-order section by adding another
///        pole and zero to the current pole and zero in the section:
///
///         1. If the current pole and zero are both complex, add their
///            conjugates.
///
///         2. Else if the pole is complex and the zero is real, add the
///            conjugate pole and the next closest real zero.
///
///         3. Else if the pole is real and the zero is complex, add the
///            conjugate zero and the real pole closest to those zeros.
///
///         4. Else (we must have a real pole and real zero) add the next
///            real pole closest to the unit circle, and then add the real
///            zero closest to that pole.
///
/// .. [#] This conditional can only be met for specific odd-order inputs
///        with the ``pairing = 'keep_odd'`` or ``'minimal'`` methods.
///
/// .. versionadded:: 0.16.0
///
/// Examples
/// --------
///
/// Design a 6th order low-pass elliptic digital filter for a system with a
/// sampling rate of 8000 Hz that has a pass-band corner frequency of
/// 1000 Hz. The ripple in the pass-band should not exceed 0.087 dB, and
/// the attenuation in the stop-band should be at least 90 dB.
///
/// In the following call to `ellip`, we could use ``output='sos'``,
/// but for this example, we'll use ``output='zpk'``, and then convert
/// to SOS format with `zpk2sos`:
///
/// >>> from scipy import signal
/// >>> z, p, k = signal.ellip(6, 0.087, 90, 1000/(0.5*8000), output='zpk')
///
/// Now convert to SOS format.
///
/// >>> sos = signal.zpk2sos(z, p, k)
///
/// The coefficients of the numerators of the sections:
///
/// >>> sos[:, :3]
/// array([[ 0.0014154 ,  0.00248707,  0.0014154 ],
///        [ 1.        ,  0.72965193,  1.        ],
///        [ 1.        ,  0.17594966,  1.        ]])
///
/// The symmetry in the coefficients occurs because all the zeros are on the
/// unit circle.
///
/// The coefficients of the denominators of the sections:
///
/// >>> sos[:, 3:]
/// array([[ 1.        , -1.32543251,  0.46989499],
///        [ 1.        , -1.26117915,  0.6262586 ],
///        [ 1.        , -1.25707217,  0.86199667]])
///
/// The next example shows the effect of the `pairing` option.  We have a
/// system with three poles and three zeros, so the SOS array will have
/// shape (2, 6). The means there is, in effect, an extra pole and an extra
/// zero at the origin in the SOS representation.
///
/// >>> z1 = np.array([-1, -0.5-0.5j, -0.5+0.5j])
/// >>> p1 = np.array([0.75, 0.8+0.1j, 0.8-0.1j])
///
/// With ``pairing='nearest'`` (the default), we obtain
///
/// >>> signal.zpk2sos(z1, p1, 1)
/// array([[ 1.  ,  1.  ,  0.5 ,  1.  , -0.75,  0.  ],
///        [ 1.  ,  1.  ,  0.  ,  1.  , -1.6 ,  0.65]])
///
/// The first section has the zeros {-0.5-0.05j, -0.5+0.5j} and the poles
/// {0, 0.75}, and the second section has the zeros {-1, 0} and poles
/// {0.8+0.1j, 0.8-0.1j}. Note that the extra pole and zero at the origin
/// have been assigned to different sections.
///
/// With ``pairing='keep_odd'``, we obtain:
///
/// >>> signal.zpk2sos(z1, p1, 1, pairing='keep_odd')
/// array([[ 1.  ,  1.  ,  0.  ,  1.  , -0.75,  0.  ],
///        [ 1.  ,  1.  ,  0.5 ,  1.  , -1.6 ,  0.65]])
///
/// The extra pole and zero at the origin are in the same section.
/// The first section is, in effect, a first-order section.
///
/// With ``pairing='minimal'``, the first-order section doesn't have
/// the extra pole and zero at the origin:
///
/// >>> signal.zpk2sos(z1, p1, 1, pairing='minimal')
/// array([[ 0.  ,  1.  ,  1.  ,  0.  ,  1.  , -0.75],
///        [ 1.  ,  1.  ,  0.5 ,  1.  , -1.6 ,  0.65]])
///
/// """
///
pub fn zpk2sos<F, const N: usize>(zpk: Zpk<F, N>, analog: bool) -> [Sos<F>; N / 2 - 1]
where
    F: RealField + Float,
    [Sos<F>; N / 2 - 1]: Sized,
{
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
