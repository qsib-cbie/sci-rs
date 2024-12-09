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
#[cfg(feature = "alloc")]
pub fn lp2bp_zpk_dyn<F>(zpk: ZpkFormatFilter<F>, wo: Option<F>, bw: Option<F>) -> ZpkFormatFilter<F>
where
    F: RealField + Float,
{
    let wo = wo.unwrap_or_else(F::one);
    let bw = bw.unwrap_or_else(F::one);

    let degree = relative_degree_dyn(&zpk.z, &zpk.p);

    // Scale poles and zeros to desired bandwidth
    let two = unsafe { F::from(2.).unwrap_unchecked() };
    let z_lp: Vec<_> = zpk
        .z
        .iter()
        .map(|zi| *zi * Complex::new(bw / two, F::zero()))
        .collect();

    let p_lp: Vec<_> = zpk
        .p
        .iter()
        .map(|zi| *zi * Complex::new(bw / two, F::zero()))
        .collect();

    // Duplicate poles and zeros and shift from baseband to +wo and -wo
    let wo2 = Complex::new(Float::powi(wo, 2), F::zero());

    let z_bp_t = z_lp
        .iter()
        .map(|zi| (*zi, (zi.powi(2) - wo2).sqrt()))
        .collect::<Vec<(Complex<F>, Complex<F>)>>();
    let mut z_bp = z_bp_t
        .iter()
        .map(|(a, b)| a + b)
        .chain(z_bp_t.iter().map(|(a, b)| a - b))
        .collect::<Vec<Complex<F>>>();

    let p_bp_t = p_lp
        .iter()
        .map(|zi| (*zi, (zi.powi(2) - wo2).sqrt()))
        .collect::<Vec<(Complex<F>, Complex<F>)>>();
    let p_bp = p_bp_t
        .iter()
        .map(|(a, b)| a + b)
        .chain(p_bp_t.iter().map(|(a, b)| a - b))
        .collect::<Vec<Complex<F>>>();

    // Move degree zeros to origin, leaving degree zeros at infinity for BPF
    z_bp.extend((0..degree).map(|_| Complex::new(F::zero(), F::zero())));

    // Cancel out gain change from frequency scaling
    let k_bp = zpk.k * Float::powi(bw, degree as i32);

    ZpkFormatFilter::new(z_bp, p_bp, k_bp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[cfg(all(feature = "alloc", feature = "std"))]
    #[test]
    fn matches_scipy_example() {
        // butter(4, [10, 50], btype='bandpass', output='sos', fs=1666)
        let zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            Vec::new(),
            vec![
                Complex::new(-0.38268343, 0.92387953),
                Complex::new(-0.92387953, 0.38268343),
                Complex::new(-0.92387953, -0.38268343),
                Complex::new(-0.38268343, -0.92387953),
            ],
            1.,
        );
        let wo = 0.16892363165758506;
        let bw = 0.30282619318434084;

        let expected_z = vec![Complex::new(0., 0.); 4];
        let expected_p = vec![
            Complex::new(-0.02022036, -0.07498294),
            Complex::new(-0.07648538, -0.06990013),
            Complex::new(-0.07648538, 0.06990013),
            Complex::new(-0.02022036, 0.07498294),
            Complex::new(-0.0956662, 0.35475786),
            Complex::new(-0.20328954, 0.1857867),
            Complex::new(-0.20328954, -0.1857867),
            Complex::new(-0.0956662, -0.35475786),
        ];
        let expected_k = 0.008409569194994788;
        let actual_zpk: ZpkFormatFilter<f64> = lp2bp_zpk_dyn(zpk, Some(wo), Some(bw));

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
