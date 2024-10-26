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

#[cfg(feature = "alloc")]
///Transform a lowpass filter prototype to a highpass filter.
///
///Return an analog high-pass filter with cutoff frequency `wo`
///from an analog low-pass filter prototype with unity cutoff frequency, in
///transfer function ('ba') representation.
///
///Parameters
///----------
///b : array_like
///   Numerator polynomial coefficients.
///a : array_like
///    Denominator polynomial coefficients.
///wo : float
///    Desired cutoff, as angular frequency (e.g., rad/s).
///    Defaults to no change.
///
///Returns
///    -------
///    z : ndarray
///        Zeros of the transformed high-pass filter transfer function.
///    p : ndarray
///        Poles of the transformed high-pass filter transfer function.
///    k : float
///        System gain of the transformed high-pass filter.
///
///See Also
///--------
///lp2lp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
///lp2hp
///
///Notes
///-----
///This is derived from the s-plane substitution
///
///.. math:: s \rightarrow \frac{\omega_0}{s}
///
///This maintains symmetry of the lowpass and highpass responses on a
///logarithmic scale.
///
///Examples
///--------
///>>> from scipy import signal
///>>> import matplotlib.pyplot as plt
///
///>>> lp = signal.lti([1.0], [1.0, 1.0])
///>>> hp = signal.lti(*signal.lp2hp(lp.num, lp.den))
///>>> w, mag_lp, p_lp = lp.bode()
///>>> w, mag_hp, p_hp = hp.bode(w)
///
///>>> plt.plot(w, mag_lp, label='Lowpass')
///>>> plt.plot(w, mag_hp, label='Highpass')
///>>> plt.semilogx()
///>>> plt.grid(True)
///>>> plt.xlabel('Frequency [rad/s]')
///>>> plt.ylabel('Magnitude \[dB\]')
///>>> plt.legend()
pub fn lp2hp_zpk_dyn<F>(zpk: ZpkFormatFilter<F>, wo: Option<F>) -> ZpkFormatFilter<F>
where
    F: RealField + Float,
{
    let wo = wo.unwrap_or_else(F::one); // Avoid int wraparound

    let degree = relative_degree_dyn(&zpk.z, &zpk.p);

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

    // # If lowpass had zeros at infinity, inverting moves them to origin.
    // z_hp = append(z_hp, zeros(degree))

    z_hp.extend((0..degree).map(|_| Complex::new(F::zero(), F::zero())));

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
        .map(|pi| Complex::new(one_neg, F::zero()) * *pi)
        .fold(Complex::new(F::one(), F::zero()), |acc, pi| acc * pi);
    let k_hp = zpk.k * (num / denom).real();

    // return z_hp, p_hp, k_hp
    ZpkFormatFilter::new(z_hp, p_hp, k_hp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[cfg(all(feature = "alloc", feature = "std"))]
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
        let expected_z = vec![Complex::new(0.16892363, 0.16892363); 4];
        let expected_p = vec![
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

    #[cfg(all(feature = "alloc", feature = "std"))]
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
        let expected_z = vec![Complex::new(0.1216, 0.1216); 4];
        let expected_p = vec![
            Complex::new(0.1307223, 0.03503196),
            Complex::new(0.15726189, 0.01745671),
            Complex::new(0.15726189, -0.01745671),
            Complex::new(0.13072233, -0.03503196),
        ];
        let expected_k = 1.4481908047355312;
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

        assert_relative_eq!(actual_zpk.k, expected_k, max_relative = 1e-6);
    }

    #[cfg(all(feature = "alloc", feature = "std"))]
    #[test]
    fn matches_scipy_example_highpass_three() {
        // >>> butter(5,[2],btype='hp',output='zpk',fs=600)
        // [] [-0.30901699+0.95105652j -0.80901699+0.58778525j -1.        -0.j
        //  -0.80901699-0.58778525j -0.30901699-0.95105652j] 1 [0.04188943]
        // [0. 0. 0. 0. 0.] [-0.01294455-0.03983922j -0.03388926-0.02462199j -0.04188943-0.j
        //  -0.03388926+0.02462199j -0.01294455+0.03983922j] 1.0
        // (array([1., 1., 1., 1., 1.]), array([0.99335214-0.01978936j, 0.98312384-0.01210456j,
        //        0.97927235+0.j        , 0.98312384+0.01210456j,
        //        0.99335214+0.01978936j]), 0.9666790028093929)

        let input_zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![],
            vec![
                Complex::new(-0.30901699, 0.95105652),
                Complex::new(-0.80901699, 0.58778525),
                Complex::new(-1., 0.),
                Complex::new(-0.80901699, -0.58778525),
                Complex::new(-0.30901699, -0.95105652),
            ],
            1.,
        );
        let wo = Some(0.04188943);

        let expected_zpk: ZpkFormatFilter<_> = ZpkFormatFilter::new(
            vec![
                Complex::new(0., 0.),
                Complex::new(0., 0.),
                Complex::new(0., 0.),
                Complex::new(0., 0.),
                Complex::new(0., 0.),
            ],
            vec![
                Complex::new(-0.01294455, -0.03983922),
                Complex::new(-0.03388926, -0.02462199),
                Complex::new(-0.04188943, 0.),
                Complex::new(-0.03388926, 0.02462199),
                Complex::new(-0.01294455, 0.03983922),
            ],
            1.0,
        );
        let actual_zpk: ZpkFormatFilter<f64> = lp2hp_zpk_dyn(input_zpk, wo);

        assert_eq!(actual_zpk.z.len(), expected_zpk.z.len());
        for (a, e) in actual_zpk.z.iter().zip(expected_zpk.z.iter()) {
            assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
            assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
        }

        assert_eq!(actual_zpk.p.len(), expected_zpk.p.len());
        for (a, e) in actual_zpk.p.iter().zip(expected_zpk.p.iter()) {
            assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
            assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
        }

        assert_relative_eq!(actual_zpk.k, expected_zpk.k, max_relative = 1e-6);
    }
}
