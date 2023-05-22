use core::{f64::consts::PI, ops::Mul};
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "use_std")]
use super::{
    relative_degree::relative_degree_dyn, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};


    // Transform a lowpass filter prototype to a different frequency.

    // Return an analog low-pass filter with cutoff frequency `wo`
    // from an analog low-pass filter prototype with unity cutoff frequency,
    // using zeros, poles, and gain ('zpk') representation.

    // Parameters
    // ----------
    // z : array_like
    //     Zeros of the analog filter transfer function.
    // p : array_like
    //     Poles of the analog filter transfer function.
    // k : float
    //     System gain of the analog filter transfer function.
    // wo : float
    //     Desired cutoff, as angular frequency (e.g., rad/s).
    //     Defaults to no change.

    // Returns
    // -------
    // z : ndarray
    //     Zeros of the transformed low-pass filter transfer function.
    // p : ndarray
    //     Poles of the transformed low-pass filter transfer function.
    // k : float
    //     System gain of the transformed low-pass filter.

    // See Also
    // --------
    // lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
    // lp2lp

    // Notes
    // -----
    // This is derived from the s-plane substitution

    // .. math:: s \rightarrow \frac{s}{\omega_0}

    // .. versionadded:: 1.1.0


#[cfg(feature = "use_std")]
pub fn lp2lp_zpk_dyn<F>(zpk: ZpkFormatFilter<F>, wo: Option<F>, bw: Option<F>) -> ZpkFormatFilter<F>
where
    F: RealField + Float,
{
    let wo = wo.unwrap_or_else(F::one);  // Avoid int wraparound

    let degree = relative_degree_dyn(&zpk.z, &zpk.p);

    // Scale all points radially from origin to shift cutoff frequency
    z_lp = wo * z;
    p_lp = wo * p;

    // Each shifted pole decreases gain by wo, each shifted zero increases it.
    // Cancel out the net change to keep overall gain the same
    k_lp = k * wo**degree;

    ZpkFormatFilter::new(z_bp, p_bp, k_bp)
} fn lp2lp_zpk_dyn