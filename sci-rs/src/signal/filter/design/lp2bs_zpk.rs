use core::{f64::consts::PI, ops::Mul};
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "alloc")]
use super::{
    relative_degree::relative_degree_dyn, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};

#[cfg(feature = "alloc")]
/// Transform a lowpass filter prototype to a bandstop filter.
///
/// Return an analog band-stop filter with center frequency `wo` and
/// stopband width `bw` from an analog low-pass filter prototype with unity
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
///     Desired stopband center, as angular frequency (e.g., rad/s).
///     Defaults to no change.
/// bw : float
///    Desired stopband width, as angular frequency (e.g., rad/s).
///     Defaults to 1.
///
/// Returns
/// -------
/// z : ndarray
///     Zeros of the transformed band-stop filter transfer function.
/// p : ndarray
///     Poles of the transformed band-stop filter transfer function.
/// k : float
///     System gain of the transformed band-stop filter.
///
/// See Also
/// --------
/// lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, bilinear
/// lp2bs
///
/// Notes
/// -----
/// This is derived from the s-plane substitution
///
//. .. math:: s \rightarrow \frac{s \cdot \mathrm{BW}}{s^2 + {\omega_0}^2}
///
/// This is the "wideband" transformation, producing a stopband with
/// geometric (log frequency) symmetry about `wo`.
pub fn lp2bs_zpk_dyn<F>(zpk: ZpkFormatFilter<F>, wo: Option<F>, bw: Option<F>) -> ZpkFormatFilter<F>
where
    F: RealField + Float,
{
    let wo = wo.unwrap_or_else(F::one);
    let bw = bw.unwrap_or_else(F::one);

    let degree = relative_degree_dyn(&zpk.z, &zpk.p);

    //Invert to a highpass filter with desired bandwidth
    let one_neg = unsafe { F::from(-1.).unwrap_unchecked() };
    let two = unsafe { F::from(2.).unwrap_unchecked() };
    let two_neg = unsafe { F::from(-2).unwrap_unchecked() };
    let z_hp: Vec<_> = zpk
        .z
        .iter()
        .map(|zi| (Complex::new(bw * one_neg, F::zero()) / *zi) * Complex::new(F::zero(), one_neg))
        .collect();

    let p_hp: Vec<_> = zpk
        .p
        .iter()
        .map(|zi| (Complex::new(F::zero(), bw / two) / *zi) * Complex::new(F::zero(), one_neg))
        .collect();

    //Duplicate poles and zeros and shift from baseband to +wo and -wo
    let wo2 = Complex::new(Float::powi(wo, 2), F::zero());

    let mut z_bs_t = z_hp
        .iter()
        .map(|zi| {
            (
                *zi,
                ((zi - Complex::new(F::zero(), zi.im)).powi(2) - wo2).sqrt(),
            )
        })
        .collect::<Vec<(Complex<F>, Complex<F>)>>();

    let z_bs_t = z_bs_t
        .iter()
        .map(|(zi, zi2)| (Complex::new(zi.re, zi2.im), Complex::new(zi.re, zi2.im)))
        .collect::<Vec<(Complex<F>, Complex<F>)>>();

    let mut z_bs = z_bs_t
        .iter()
        .map(|(a, b)| Complex::new(a.re, b.im))
        .chain(
            z_bs_t
                .iter()
                .map(|(a, b)| Complex::new(a.re, b.im * one_neg)),
        )
        .collect::<Vec<Complex<F>>>();

    let p_bs_t = p_hp
        .iter()
        .map(|zi| (*zi, (zi.powi(2) - wo2).sqrt()))
        .collect::<Vec<(Complex<F>, Complex<F>)>>();

    let p_bs = p_bs_t
        .iter()
        .map(|(a, b)| a + b)
        .chain(p_bs_t.iter().map(|(a, b)| a - b))
        .collect::<Vec<Complex<F>>>();

    //Move any zeros that were at infinity to the center of the stopband
    z_bs.extend((0..degree).map(|_| Complex::new(F::zero(), F::one() * wo)));
    z_bs.extend((0..degree).map(|_| Complex::new(F::zero(), one_neg * wo)));

    //Cancel out gain change caused by inversion
    let num = zpk.z.iter().copied().fold(F::one(), |acc, zi| acc * zi.re);
    let denom = zpk
        .p
        .iter()
        .map(|pi| Complex::new(one_neg, F::zero()) * *pi)
        .fold(Complex::new(F::one(), F::zero()), |acc, pi| acc * pi);
    let k_bs = zpk.k * (num / denom.real());

    ZpkFormatFilter::new(z_bs, p_bs, k_bs)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[cfg(feature = "alloc")]
    #[test]
    fn matches_scipy_example_bandstop() {
        let zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![
                Complex::new(-1., -1.),
                Complex::new(-1., -1.),
                Complex::new(-1., -1.),
                Complex::new(-1., -1.),
            ],
            vec![
                Complex::new(0.98765384, 0.02863265),
                Complex::new(0.97136157, 0.01166439),
                Complex::new(0.97136157, -0.01166439),
                Complex::new(0.98765384, -0.02863265),
            ],
            5.8105542410017214e-08,
        );

        let wo = 1.6;
        let bw = 1.2;
        let expected_z: Vec<_> = vec![
            Complex::new(-0.6, -1.4832397),
            Complex::new(-0.6, -1.4832397),
            Complex::new(-0.6, -1.4832397),
            Complex::new(-0.6, -1.4832397),
            Complex::new(-0.6, 1.4832397),
            Complex::new(-0.6, 1.4832397),
            Complex::new(-0.6, 1.4832397),
            Complex::new(-0.6, 1.4832397),
        ];
        let expected_p: Vec<_> = vec![
            Complex::new(0.61420466, -1.49811199),
            Complex::new(0.62070377, -1.48343601),
            Complex::new(0.62070377, 1.48343601),
            Complex::new(0.61420466, 1.49811199),
            Complex::new(0.59977563, 1.46291801),
            Complex::new(0.61449745, 1.46860336),
            Complex::new(0.61449745, -1.46860336),
            Complex::new(0.59977563, -1.46291801),
        ];
        let expected_k = 6.306940631261856e-08;
        let actual_zpk: ZpkFormatFilter<f64> = lp2bs_zpk_dyn(zpk, Some(wo), Some(bw));

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

    #[cfg(feature = "alloc")]
    #[test]
    fn matches_scipy_example_bandstop_two() {
        let zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![
                Complex::new(1., 1.),
                Complex::new(1., 1.),
                Complex::new(1., 1.),
                Complex::new(1., 1.),
            ],
            vec![
                Complex::new(0.86788666, -0.23258286),
                Complex::new(0.76382075, -0.08478723),
                Complex::new(0.76382075, 0.08478723),
                Complex::new(0.86788666, 0.23258286),
            ],
            0.6905166297398233,
        );

        let wo = 1.6;
        let bw = 1.2;
        let expected_z: Vec<_> = vec![
            Complex::new(0.6, 1.4832397),
            Complex::new(0.6, 1.4832397),
            Complex::new(0.6, 1.4832397),
            Complex::new(0.6, 1.4832397),
            Complex::new(0.6, -1.4832397),
            Complex::new(0.6, -1.4832397),
            Complex::new(0.6, -1.4832397),
            Complex::new(0.6, -1.4832397),
        ];
        let expected_p: Vec<_> = vec![
            Complex::new(0.72053235, 1.64918244),
            Complex::new(0.82361252, 1.48883633),
            Complex::new(0.82361252, -1.48883633),
            Complex::new(0.72053235, -1.64918244),
            Complex::new(0.56949063, -1.30347227),
            Complex::new(0.728314, -1.31656612),
            Complex::new(0.728314, 1.31656612),
            Complex::new(0.56949063, 1.30347227),
        ];
        let expected_k = 1.4481908047355312;
        let actual_zpk: ZpkFormatFilter<f64> = lp2bs_zpk_dyn(zpk, Some(wo), Some(bw));

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

        assert_relative_eq!(actual_zpk.k, expected_k, max_relative = 1e-6);
    }
}
