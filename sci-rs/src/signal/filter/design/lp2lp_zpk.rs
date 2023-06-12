use core::{f64::consts::PI, ops::Mul};
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "use_std")]
use super::{
    relative_degree::relative_degree_dyn, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};


    /// Transform a lowpass filter prototype to a different frequency.

    /// Return an analog low-pass filter with cutoff frequency `wo`
    /// from an analog low-pass filter prototype with unity cutoff frequency,
    /// using zeros, poles, and gain ('zpk') representation.

    /// Parameters
    /// ----------
    /// z : array_like
    ///     Zeros of the analog filter transfer function.
    /// p : array_like
    ///     Poles of the analog filter transfer function.
    /// k : float
    ///     System gain of the analog filter transfer function.
    /// wo : float
    ///     Desired cutoff, as angular frequency (e.g., rad/s).
    ///     Defaults to no change.

    /// Returns
    /// -------
    /// z : ndarray
    ///     Zeros of the transformed low-pass filter transfer function.
    /// p : ndarray
    ///     Poles of the transformed low-pass filter transfer function.
    /// k : float
    ///     System gain of the transformed low-pass filter.

    /// See Also
    /// --------
    /// lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
    /// lp2lp

    /// Notes
    /// -----
    /// This is derived from the s-plane substitution

    /// .. math:: s \rightarrow \frac{s}{\omega_0}

    /// .. versionadded:: 1.1.0


#[cfg(feature = "use_std")]
pub fn lp2lp_zpk_dyn<F>(zpk: ZpkFormatFilter<F>, wo: Option<F>) -> ZpkFormatFilter<F>
where
    F: RealField + Float,
{
    let wo = wo.unwrap_or_else(F::one); // Avoid int wraparound

    let degree = relative_degree_dyn(&zpk.z, &zpk.p);

    /// Scale all points radially from origin to shift cutoff frequency
    let z_lp: Vec<_> = zpk
        .z
        .iter()
        .map(|zi| *zi * Complex::new(wo, F::zero()))
        .collect();

    let p_lp: Vec<_> = zpk
        .p
        .iter()
        .map(|zi| *zi * Complex::new(wo, F::zero()))
        .collect();


    /// Each shifted pole decreases gain by wo, each shifted zero increases it.
    /// Cancel out the net change to keep overall gain the same

    let k_lp = zpk.k * Float::powi(wo, degree as i32);

    ZpkFormatFilter::new(z_lp, p_lp, k_lp)
} fn lp2lp_zpk_dyn