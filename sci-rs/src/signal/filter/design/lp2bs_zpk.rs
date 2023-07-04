use core::{f64::consts::PI, ops::Mul};
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "use_std")]
use super::{
    relative_degree::relative_degree_dyn, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};

#[cfg(feature = "use_std")]
///Transform a lowpass filter prototype to a bandstop filter.
///
///Return an analog band-stop filter with center frequency `wo` and
///stopband width `bw` from an analog low-pass filter prototype with unity
///cutoff frequency, using zeros, poles, and gain ('zpk') representation.
///
///Parameters
///----------
///z : array_like
///    Zeros of the analog filter transfer function.
///p : array_like
///    Poles of the analog filter transfer function.
///k : float
///    System gain of the analog filter transfer function.
///wo : float
///    Desired stopband center, as angular frequency (e.g., rad/s).
///    Defaults to no change.
///bw : float
///   Desired stopband width, as angular frequency (e.g., rad/s).
///    Defaults to 1.
///
///Returns
///-------
///z : ndarray
///    Zeros of the transformed band-stop filter transfer function.
///p : ndarray
///    Poles of the transformed band-stop filter transfer function.
///k : float
///    System gain of the transformed band-stop filter.
///
///See Also
///--------
///lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, bilinear
///lp2bs
///
///Notes
///-----
///This is derived from the s-plane substitution
///
//... math:: s \rightarrow \frac{s \cdot \mathrm{BW}}{s^2 + {\omega_0}^2}
///
///This is the "wideband" transformation, producing a stopband with
///geometric (log frequency) symmetry about `wo`.
pub fn lp2bs_zpk_dyn<F>(zpk: ZpkFormatFilter<F>, wo: Option<F>, bw: Option<F>) -> ZpkFormatFilter<F>
where
    F: RealField + Float,
{
    // wo = float(wo)
    // bw = float(bw)

    let wo = wo.unwrap_or_else(F::one); // Avoid int wraparound

    // degree = _relative_degree(z, p)

    let degree = relative_degree_dyn(&zpk.z, &zpk.p);

    // # Invert to a highpass filter with desired bandwidth
    // z_hp = (bw/2) / z
    // p_hp = (bw/2) / p

    // # Invert positions radially about unit circle to convert LPF to HPF
    // # Scale all points radially from origin to shift cutoff frequency
    // z_hp = wo / z
    // p_hp = wo / p
    let one_neg = unsafe { F::from(-1.).unwrap_unchecked() };
    let two = unsafe { F::from(-2.).unwrap_unchecked() };
    let mut z_hp: Vec<_> = zpk
        .z
        .iter()
        .map(|zi| (Complex::new(wo * two, F::zero()) / *zi) * (Complex::new(F::zero(), one_neg)))
        .collect();

    let p_hp: Vec<_> = zpk
        .p
        .iter()
        .map(|zi| Complex::new(wo, F::zero()) / *zi)
        .collect();

    // # Square root needs to produce complex result, not NaN
    // z_hp = z_hp.astype(complex)
    // p_hp = p_hp.astype(complex)

    // # Duplicate poles and zeros and shift from baseband to +wo and -wo
    // z_bs = concatenate((z_hp + sqrt(z_hp**2 - wo**2),
    //                     z_hp - sqrt(z_hp**2 - wo**2)))
    // p_bs = concatenate((p_hp + sqrt(p_hp**2 - wo**2),
    //                     p_hp - sqrt(p_hp**2 - wo**2)))

    // # Move any zeros that were at infinity to the center of the stopband
    // z_bs = append(z_bs, full(degree, +1j*wo))
    // z_bs = append(z_bs, full(degree, -1j*wo))

    // # If lowpass had zeros at infinity, inverting moves them to origin.
    // z_hp = append(z_hp, zeros(degree))

    z_hp.extend((0..degree).map(|_| Complex::new(F::zero(), F::zero())));

    // # Cancel out gain change caused by inversion
    // k_bs = k * real(prod(-z) / prod(-p))

    // # Cancel out gain change caused by inversion
    // k_hp = k * real(prod(-z) / prod(-p))

    let num = zpk
        .z
        .iter()
        .map(|zi| Complex::new(one_neg, F::zero()) + *zi)
        .fold(Complex::new(F::one(), F::zero()), |acc, zi| acc * zi);
    let denom = zpk
        .p
        .iter()
        .map(|pi| Complex::new(F::one(), F::zero()) * *pi)
        .fold(Complex::new(F::one(), F::zero()), |acc, pi| acc * pi);
    let k_hp = zpk.k * (num / denom).real();

    // return z_bs, p_bs, k_bs

    // return z_hp, p_hp, k_hp
    ZpkFormatFilter::new(z_hp, p_hp, k_hp)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_example_highpass() {
        //zo = [1. 1. 1. 1.]
        //po = [0.91652519-0.16159409j 0.83726957-0.06114637j 0.83726957+0.06114637j 0.91652519+0.16159409j]
        //ko = 0.7812898609101183

        let zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![
                Complex::new(1., 1.),
                Complex::new(1., 1.),
                Complex::new(1., 1.),
                Complex::new(1., 1.),
            ],
            vec![
                Complex::new(0.91652519, -0.16159409),
                Complex::new(0.83726957, -0.06114637),
                Complex::new(0.83726957, 0.06114637),
                Complex::new(0.91652519, 0.16159409),
            ],
            0.7812898609101183,
        );

        let wo = 0.16892363165758506;
        //z1 = [0.16892363 0.16892363 0.16892363 0.16892363]
        //p1 = [0.17875212+0.03151608j 0.20068502+0.20068502j 0.20068502-0.01465616j 0.17875212-0.03151608j]
        //k1 = 1.2799346870309942
        let expected_z: Vec<_> = vec![Complex::new(0.16892363, 0.16892363); 4];
        let expected_p: Vec<_> = vec![
            Complex::new(0.17875212, 0.03151608),
            Complex::new(0.20068502, 0.01465616),
            Complex::new(0.20068502, -0.01465616),
            Complex::new(0.17875212, -0.03151608),
        ];
        let expected_k = 1.2799346870309942;
        let actual_zpk: ZpkFormatFilter<f64> = lp2hp_zpk_dyn(zpk, Some(wo));

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
    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_example_highpass_two() {
        //zo = [1. 1. 1. 1.]
        //po = [0.86788666-0.23258286j 0.76382075-0.08478723j 0.76382075+0.08478723j 0.86788666+0.23258286j]
        //ko = 0.6905166297398233

        let zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![
                Complex::new(1., 1.),
                Complex::new(1., 1.),
                Complex::new(1., 1.),
                Complex::new(1., 1.),
            ],
            vec![
                Complex::new(0.86788666, -0.23258286),
                Complex::new(0.76382075, -0.084787237),
                Complex::new(0.76382075, 0.08478723),
                Complex::new(0.86788666, 0.23258286),
            ],
            0.6905166297398233,
        );

        let wo = 0.1216;
        //z1 = [0.1216 0.1216 0.1216 0.1216]
        //p1 = [0.13072233+0.03503196j 0.15726189+0.01745671j 0.15726189-0.01745671j 0.13072233-0.03503196j]
        //k1 = 1.4481908047355312
        let expected_z: Vec<_> = vec![Complex::new(0.1216, 0.1216); 4];
        let expected_p: Vec<_> = vec![
            Complex::new(0.1307223, 0.03503196),
            Complex::new(0.15726189, 0.01745671),
            Complex::new(0.15726189, -0.01745671),
            Complex::new(0.13072233, -0.03503196),
        ];
        let expected_k = 1.4481908047355312;
        let actual_zpk: ZpkFormatFilter<f64> = lp2bs_zpk_dyn(zpk, Some(wo));

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
