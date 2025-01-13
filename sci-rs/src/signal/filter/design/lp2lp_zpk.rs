use core::{f64::consts::PI, ops::Mul};
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "alloc")]
use super::{
    relative_degree::relative_degree_dyn, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

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
#[cfg(feature = "alloc")]
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[cfg(all(feature = "alloc", feature = "std"))]
    #[test]
    fn matches_scipy_example_lowpass() {
        // signal.butter(4, 10, btype='lowpass', output='sos', fs=1666)
        // zo = [-1. -1. -1. -1.]
        // po = [0.98507508+0.03433962j 0.96565036+0.01394346j 0.96565036-0.01394346j  0.98507508-0.03433962j]
        //ko = 1.204213960778651e-07
        let zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![
                Complex::new(-1., -1.),
                Complex::new(-1., -1.),
                Complex::new(-1., -1.),
                Complex::new(-1., -1.),
            ],
            vec![
                Complex::new(0.98507508, 0.03433962),
                Complex::new(0.96565036, 0.01394346),
                Complex::new(0.96565036, -0.01394346),
                Complex::new(0.98507508, -0.03433962),
            ],
            1.204213960778651e-07,
        );

        let wo = 0.16892363165758506;
        // z1 = [-0.16892363 -0.16892363 -0.16892363 -0.16892363]
        // p1 = [0.16640246+0.00580077j 0.16312117+0.00235538j 0.16312117-0.00235538j 0.16640246-0.00580077j]
        // k1 = 1.204213960778651e-07
        let expected_z = vec![Complex::new(-0.16892363, -0.16892363); 4];
        let expected_p = vec![
            Complex::new(0.16640246, 0.00580077),
            Complex::new(0.16312117, 0.00235538),
            Complex::new(0.16312117, -0.00235538),
            Complex::new(0.16640246, -0.00580077),
        ];
        let expected_k = 1.204213960778651e-07;
        let actual_zpk: ZpkFormatFilter<f64> = lp2lp_zpk_dyn(zpk, Some(wo));

        // println!("BAHH: {}", actual_zpk.k);

        assert_eq!(actual_zpk.z.len(), expected_z.len());
        for (a, e) in actual_zpk.z.iter().zip(expected_z.iter()) {
            assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
            assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
        }

        assert_eq!(actual_zpk.p.len(), expected_p.len());
        for (a, e) in actual_zpk.p.iter().zip(expected_p.iter()) {
            assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
            assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
        }

        assert_relative_eq!(actual_zpk.k, expected_k, max_relative = 1e-8);
    }
}
