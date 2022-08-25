use core::{f64::consts::PI, ops::Mul};

use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{FilterBandType, FilterOutput, FilterOutputType, FilterType, Sos, Zpk};

/// """
/// Transform a lowpass filter prototype to a bandpass filter.
///
/// Return an analog band-pass filter with center frequency `wo` and
/// bandwidth `bw` from an analog low-pass filter prototype with unity
/// cutoff frequency, using zeros, poles, and gain ('zpk') representation.
///
/// Parameters
/// ----------
/// z : array_like
///     Zeros of the analog filter transfer function.
/// p : array_like
///     Poles of the analog filter transfer function.
/// k : float
///     System gain of the analog filter transfer function.
/// wo : float
///     Desired passband center, as angular frequency (e.g., rad/s).
///     Defaults to no change.
/// bw : float
///     Desired passband width, as angular frequency (e.g., rad/s).
///     Defaults to 1.
///
/// Returns
/// -------
/// z : ndarray
///     Zeros of the transformed band-pass filter transfer function.
/// p : ndarray
///     Poles of the transformed band-pass filter transfer function.
/// k : float
///     System gain of the transformed band-pass filter.
///
/// See Also
/// --------
/// lp2lp_zpk, lp2hp_zpk, lp2bs_zpk, bilinear
/// lp2bp
///
/// Notes
/// -----
/// This is derived from the s-plane substitution
///
/// .. math:: s \rightarrow \frac{s^2 + {\omega_0}^2}{s \cdot \mathrm{BW}}
///
/// This is the "wideband" transformation, producing a passband with
/// geometric (log frequency) symmetry about `wo`.
///
/// .. versionadded:: 1.1.0
///
/// """
pub fn lp2bp_zpk<F, const N: usize>(zpk: Zpk<F, N>, wo: Option<F>, bw: Option<F>) -> Zpk<F, N>
where
    F: RealField + Float,
    [Sos<F>; N / 2 - 1]: Sized,
{
    let wo = wo.unwrap_or(F::one());
    let bw = bw.unwrap_or(F::one());

    // z = atleast_1d(z)
    // p = atleast_1d(p)
    // wo = float(wo)
    // bw = float(bw)

    // degree = _relative_degree(z, p)

    // # Scale poles and zeros to desired bandwidth
    // z_lp = z * bw/2
    // p_lp = p * bw/2

    // # Square root needs to produce complex result, not NaN
    // z_lp = z_lp.astype(complex)
    // p_lp = p_lp.astype(complex)

    // # Duplicate poles and zeros and shift from baseband to +wo and -wo
    // z_bp = concatenate((z_lp + sqrt(z_lp**2 - wo**2),
    //                     z_lp - sqrt(z_lp**2 - wo**2)))
    // p_bp = concatenate((p_lp + sqrt(p_lp**2 - wo**2),
    //                     p_lp - sqrt(p_lp**2 - wo**2)))

    // # Move degree zeros to origin, leaving degree zeros at infinity for BPF
    // z_bp = append(z_bp, zeros(degree))

    // # Cancel out gain change from frequency scaling
    // k_bp = k * bw**degree

    // return z_bp, p_bp, k_bp
    todo!()
}
