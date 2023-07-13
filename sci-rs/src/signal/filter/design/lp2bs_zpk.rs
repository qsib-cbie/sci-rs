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
    let wo = wo.unwrap_or_else(F::one); // Avoid int wraparound
    let bw = bw.unwrap_or_else(F::one); // Avoid int wraparound

    // degree = _relative_degree(z, p)

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

    println!("z_hp = {:?}", z_hp);
    println!("p_hp = {:?}", p_hp);
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

    println!("z_bs_t = {:?}", z_bs_t);

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

    println!("p_bs_t = {:?}", p_bs_t);

    let p_bs = p_bs_t
        .iter()
        .map(|(a, b)| a + b)
        .chain(p_bs_t.iter().map(|(a, b)| a - b))
        .collect::<Vec<Complex<F>>>();

    println!("z_bs = {:?}", z_bs);
    println!("p_bs = {:?}", p_bs);
    //Move any zeros that were at infinity to the center of the stopband
    z_bs.extend((0..degree).map(|_| Complex::new(F::zero(), F::one() * wo)));
    println!("z_bs = {:?}", z_bs);
    z_bs.extend((0..degree).map(|_| Complex::new(F::zero(), one_neg * wo)));
    println!("z_bs = {:?}", z_bs);

    //Cancel out gain change caused by inversion
    println!("z = {:?}", zpk.z);
    let num = zpk.z.iter().copied().fold(F::one(), |acc, zi| acc * zi.re);
    let denom = zpk
        .p
        .iter()
        .map(|pi| Complex::new(one_neg, F::zero()) * *pi)
        .fold(Complex::new(F::one(), F::zero()), |acc, pi| acc * pi);
    println!("num = {:?}", num);
    println!("denom = {:?}", denom);
    let k_bs = zpk.k * (num / denom.real());

    // return z_bs, p_bs, k_bs
    ZpkFormatFilter::new(z_bs, p_bs, k_bs)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_example_bandstop() {
        // zo = [-1. -1. -1. -1.]
        //po = [0.98765384+0.02863265j 0.97136157+0.01166439j 0.97136157-0.01166439j 0.98765384-0.02863265j]
        //ko = 5.8105542410017214e-08
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
        // z1 = [-0.6-1.4832397j -0.6-1.4832397j -0.6-1.4832397j -0.6-1.4832397j -0.6+1.4832397j -0.6+1.4832397j -0.6+1.4832397j -0.6+1.4832397j]
        //p1 = [0.61420466-1.49811199j 0.62070377-1.48343601j 0.62070377+1.48343601j 0.61420466+1.49811199j 0.59977563+1.46291801j 0.61449745+1.46860336j 0.61449745-1.46860336j 0.59977563-1.46291801j]
        //k1 = 6.306940631261856e-08
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
    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_example_bandstop_two() {
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
                Complex::new(0.76382075, -0.08478723),
                Complex::new(0.76382075, 0.08478723),
                Complex::new(0.86788666, 0.23258286),
            ],
            0.6905166297398233,
        );

        let wo = 1.6;
        let bw = 1.2;
        // z1 = [0.6+1.4832397j 0.6+1.4832397j 0.6+1.4832397j 0.6+1.4832397j 0.6-1.4832397j 0.6-1.4832397j 0.6-1.4832397j 0.6-1.4832397j]
        //p1 = [0.72053235+1.64918244j 0.82361252+1.48883633j 0.82361252-1.48883633j 0.72053235-1.64918244j 0.56949063-1.30347227j 0.728314  -1.31656612j 0.728314  +1.31656612j 0.56949063+1.30347227j]
        //k1 = 1.4481908047355312
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
