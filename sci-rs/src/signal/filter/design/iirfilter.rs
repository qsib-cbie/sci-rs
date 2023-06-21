use core::{f64::consts::PI, iter::Sum, ops::Mul};

#[cfg(feature = "unstable")]
use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

#[cfg(feature = "use_std")]
use crate::signal::filter::design::{zpk2tf_dyn, ZpkFormatFilter};
#[cfg(feature = "unstable")]
use crate::signal::filter::design::{zpk2tf_st, ZpkFormatFilter};

#[cfg(feature = "use_std")]
use super::{
    bilinear_zpk_dyn, lp2bp_zpk_dyn, lp2lp_zpk_dyn, zpk2sos_dyn, DigitalFilter, FilterBandType,
    FilterOutputType, FilterType, Sos,
};

///
///
/// IIR digital and analog filter design given order and critical points.
///
/// Design an Nth-order digital or analog filter and return the filter
/// coefficients.
///
/// -------
/// b, a : ndarray, ndarray
///     Numerator ('b') and denominator ('a') polynomials of the IIR filter.
///     Only returned if 'output='ba'.
/// z, p, k : ndarray, ndarray, float
///     Zeros, poles, and system gain of the IIR filter transfer
///     function.  Only returned if 'output='zpk'.
/// sos : ndarray
///     Second-order sections representation of the IIR filter.
///     Only returned if 'output=='sos'.
///
/// See Also
/// --------
/// butter : Filter design using order and critical points
/// cheby1, cheby2, ellip, bessel
/// buttord : Find order and critical points from passband and stopband spec
/// cheb1ord, cheb2ord, ellipord
/// iirdesign : General filter design using passband and stopband spec
///
#[allow(clippy::too_many_arguments)]
#[cfg(feature = "unstable")]
pub fn iirfilter_st<F, const N: usize, const M: usize>(
    wn: [F; M],
    rp: Option<F>,
    rs: Option<F>,
    btype: Option<FilterBandType>,
    ftype: Option<FilterType>,
    analog: Option<bool>,
    output: Option<FilterOutputType>,
    fs: Option<F>,
) -> DigitalFilter<F, N>
where
    F: RealField + Float + Sum,
    [(); { N * 2 + 1 }]: Sized,
    [(); { N * 2 }]: Sized,
{
    let analog = analog.unwrap_or(false);
    let mut wn = wn;

    if let Some(fs) = fs {
        if analog {
            panic!("fs cannot be specified for an analog filter");
        }

        wn.iter_mut().for_each(|wni| {
            *wni = F::from(2.).unwrap() * *wni / fs;
        });
    }

    if wn.iter().any(|wi| *wi <= F::zero()) {
        panic!("filter critical frequencies must be greater than 0");
    }

    if M > 1 && wn[0] >= wn[1] {
        panic!("Wn[0] must be less than Wn[1]");
    }

    if let Some(rp) = rp {
        if rp < F::zero() {
            panic!("passband ripple (rp) must be positive");
        }
    }

    if let Some(rs) = rs {
        if rs < F::zero() {
            panic!("stopband attenuation (rs) must be positive");
        }
    }

    // Get analog lowpass prototype
    let ftype = ftype.unwrap_or(FilterType::Butterworth);
    let zpk: ZpkFormatFilter<F, N> = match ftype {
        FilterType::Butterworth => buttap_st::<F, N>(),
        FilterType::ChebyshevI => {
            if rp.is_none() {
                panic!("passband ripple (rp) must be provided to design a Chebyshev I filter");
            }
            // cheb1ap::<N>(rp)
            todo!()
        }
        FilterType::ChebyshevII => {
            if rp.is_none() {
                panic!(
                    "stopband attenuation (rs) must be provided to design an Chebyshev II filter."
                );
            }
            // cheb2ap::<N>(rs)
            todo!()
        }
        FilterType::CauerElliptic => {
            if rs.is_none() || rp.is_none() {
                panic!("Both rp and rs must be provided to design an elliptic filter.");
            }
            // ellipap::<N>(rp, rs)
            todo!()
        }
        FilterType::BesselThomson(norm) => {
            // besselap::<N>(norm = norm),
            todo!()
        }
    };

    // Pre-warp frequencies for digital filter design
    let (fs, warped) = if !analog {
        if wn.iter().any(|wi| *wi <= F::zero() || *wi >= F::one()) {
            if let Some(fs) = fs {
                panic!(
                    "Digital filter critical frequencies must be 0 < Wn < fs/2 (fs={} -> fs/2={})",
                    fs,
                    fs / F::from(2.).unwrap()
                );
            }
            panic!("Digital filter critical frequencies must be 0 < Wn < 1");
        }
        let fs = F::from(2.).unwrap();
        let mut warped = [F::zero(); M];
        for i in 0..M {
            warped[i] = F::from(2.).unwrap() * fs * Float::tan(F::from(PI).unwrap() * wn[i] / fs);
        }
        (fs, warped)
    } else {
        (fs.unwrap_or_else(F::one), wn)
    };

    // transform to lowpass, bandpass, highpass, or bandstop
    let btype = btype.unwrap_or(FilterBandType::Bandpass);
    let zpk = match btype {
        FilterBandType::Lowpass => {
            if M != 1 {
                panic!(
                    "Must specify a single critical frequency Wn for lowpass or highpass filter"
                );
            }

            // lp2lp_zpk(z, p, k, wo = warped)
            todo!()
        }
        FilterBandType::Highpass => {
            if M != 1 {
                panic!(
                    "Must specify a single critical frequency Wn for lowpass or highpass filter"
                );
            }

            // lp2hp_zpk(z, p, k, wo = warped)
            todo!()
        }
        FilterBandType::Bandpass => {
            if M != 2 {
                panic!(
                    "Wn must specify start and stop frequencies for bandpass or bandstop filter"
                );
            }

            let bw = warped[1] - warped[0];
            let wo = Float::sqrt(warped[0] * warped[1]);
            lp2bp_zpk_st::<F, N, { N * 2 }>(zpk, Some(wo), Some(bw))
        }
        FilterBandType::Bandstop => {
            if M != 2 {
                panic!(
                    "Wn must specify start and stop frequencies for bandpass or bandstop filter"
                );
            }

            let bw = warped[1] - warped[0];
            let wo = Float::sqrt(warped[0] * warped[1]);
            // lp2bs_zpk(zpk, Some(wo), Some(bw))
            todo!()
        }
    };

    // Find discrete equivalent if necessary
    let zpk = if !analog { bilinear_zpk(zpk, fs) } else { zpk };

    // Transform to proper out type (pole-zero, state-space, numer-denom)
    let output = output.unwrap_or(FilterOutputType::Ba);
    match output {
        FilterOutputType::Zpk => DigitalFilter::Zpk(zpk),
        FilterOutputType::Ba => DigitalFilter::Ba(zpk2tf_st(&zpk.z, &zpk.p, zpk.k)),
        FilterOutputType::Sos => DigitalFilter::Sos(zpk2sos_st(zpk, None, Some(analog))),
    }
}

/// """Return (z,p,k) for analog prototype of Nth-order Butterworth filter.
///
/// The filter will have an angular (e.g., rad/s) cutoff frequency of 1.
///
/// See Also
/// --------
/// butter : Filter design function using this prototype
///
/// """
#[cfg(feature = "unstable")]
fn buttap_st<F, const N: usize>() -> ZpkFormatFilter<F, N>
where
    F: Float + RealField + Mul<Output = F>,
{
    let p: heapless::Vec<Complex<F>, N> = ((-(N as isize) + 1)..N as isize)
        .step_by(2)
        .map(|i| {
            let mi = F::from(i).unwrap();
            let num = unsafe { F::from(PI).unwrap_unchecked() * mi };
            let denom = unsafe { F::from(2. * N as f64).unwrap_unchecked() };
            let c = Complex::new(F::zero(), num / denom);
            -c.exp()
        })
        .collect::<heapless::Vec<_, N>>();

    ZpkFormatFilter::new(heapless::Vec::new(), p, F::one())
}

#[cfg(feature = "use_std")]
#[allow(clippy::too_many_arguments)]
pub fn iirfilter_dyn<F>(
    order: usize,
    wn: Vec<F>,
    rp: Option<F>,
    rs: Option<F>,
    btype: Option<FilterBandType>,
    ftype: Option<FilterType>,
    analog: Option<bool>,
    output: Option<FilterOutputType>,
    fs: Option<F>,
) -> DigitalFilter<F>
where
    F: RealField + Float + Sum,
{
    use super::bilinear_zpk_dyn;

    let analog = analog.unwrap_or(false);
    let mut wn = wn;

    if wn.len() > 2 {
        panic!("Wn may be of len 1 or 2");
    }

    if let Some(fs) = fs {
        if analog {
            panic!("fs cannot be specified for an analog filter");
        }

        wn.iter_mut().for_each(|wni| {
            *wni = F::from(2.).unwrap() * *wni / fs;
        });
    }

    if wn.iter().any(|wi| *wi <= F::zero()) {
        panic!("filter critical frequencies must be greater than 0");
    }

    if wn.len() > 1 && wn[0] >= wn[1] {
        panic!("Wn[0] must be less than Wn[1]");
    }

    if let Some(rp) = rp {
        if rp < F::zero() {
            panic!("passband ripple (rp) must be positive");
        }
    }

    if let Some(rs) = rs {
        if rs < F::zero() {
            panic!("stopband attenuation (rs) must be positive");
        }
    }

    // Get analog lowpass prototype
    let ftype = ftype.unwrap_or(FilterType::Butterworth);
    let zpk: ZpkFormatFilter<F> = match ftype {
        FilterType::Butterworth => buttap_dyn(order),
        FilterType::ChebyshevI => {
            if rp.is_none() {
                panic!("passband ripple (rp) must be provided to design a Chebyshev I filter");
            }
            // cheb1ap::<N>(rp)
            todo!()
        }
        FilterType::ChebyshevII => {
            if rp.is_none() {
                panic!(
                    "stopband attenuation (rs) must be provided to design an Chebyshev II filter."
                );
            }
            // cheb2ap::<N>(rs)
            todo!()
        }
        FilterType::CauerElliptic => {
            if rs.is_none() || rp.is_none() {
                panic!("Both rp and rs must be provided to design an elliptic filter.");
            }
            // ellipap::<N>(rp, rs)
            todo!()
        }
        FilterType::BesselThomson(norm) => {
            // besselap::<N>(norm = norm),
            todo!()
        }
    };

    // Pre-warp frequencies for digital filter design
    let (fs, warped) = if !analog {
        if wn.iter().any(|wi| *wi <= F::zero() || *wi >= F::one()) {
            if let Some(fs) = fs {
                panic!(
                    "Digital filter critical frequencies must be 0 < Wn < fs/2 (fs={} -> fs/2={})",
                    fs,
                    fs / F::from(2.).unwrap()
                );
            }
            panic!("Digital filter critical frequencies must be 0 < Wn < 1");
        }
        let fs = F::from(2.).unwrap();
        let mut warped = wn
            .iter()
            .map(|wni| F::from(2.).unwrap() * fs * Float::tan(F::from(PI).unwrap() * *wni / fs))
            .collect::<Vec<_>>();
        (fs, warped)
    } else {
        (fs.unwrap_or_else(F::one), wn.clone())
    };

    // transform to lowpass, bandpass, highpass, or bandstop
    let btype = btype.unwrap_or(FilterBandType::Bandpass);
    let zpk = match btype {
        FilterBandType::Lowpass => {
            if wn.len() != 1 {
                panic!(
                    "Must specify a single critical frequency Wn for lowpass or highpass filter"
                );
            }

            let wo = Float::sqrt(warped[0]);
            lp2lp_zpk_dyn(zpk, Some(wo))
        }
        FilterBandType::Highpass => {
            if wn.len() != 1 {
                panic!(
                    "Must specify a single critical frequency Wn for lowpass or highpass filter"
                );
            }

            todo!()
        }
        FilterBandType::Bandpass => {
            if wn.len() != 2 {
                panic!(
                    "Wn must specify start and stop frequencies for bandpass or bandstop filter"
                );
            }

            let bw = warped[1] - warped[0];
            let wo = Float::sqrt(warped[0] * warped[1]);
            lp2bp_zpk_dyn(zpk, Some(wo), Some(bw))
        }
        FilterBandType::Bandstop => {
            if wn.len() != 2 {
                panic!(
                    "Wn must specify start and stop frequencies for bandpass or bandstop filter"
                );
            }

            let bw = warped[1] - warped[0];
            let wo = Float::sqrt(warped[0] * warped[1]);
            // lp2bs_zpk(zpk, Some(wo), Some(bw))
            todo!()
        }
    };

    // Find discrete equivalent if necessary
    let zpk = if !analog {
        bilinear_zpk_dyn(zpk, fs)
    } else {
        zpk
    };

    // Transform to proper out type (pole-zero, state-space, numer-denom)
    let output = output.unwrap_or(FilterOutputType::Ba);
    match output {
        FilterOutputType::Zpk => DigitalFilter::Zpk(zpk),
        FilterOutputType::Ba => DigitalFilter::Ba(zpk2tf_dyn(2 * order, &zpk.z, &zpk.p, zpk.k)),
        FilterOutputType::Sos => DigitalFilter::Sos(zpk2sos_dyn(order, zpk, None, Some(analog))),
    }
}

/// """Return (z,p,k) for analog prototype of Nth-order Butterworth filter.
///
/// The filter will have an angular (e.g., rad/s) cutoff frequency of 1.
///
/// See Also
/// --------
/// butter : Filter design function using this prototype
///
/// """
#[cfg(feature = "use_std")]
fn buttap_dyn<F>(order: usize) -> ZpkFormatFilter<F>
where
    F: Float + RealField + Mul<Output = F>,
{
    let p: Vec<Complex<F>> = ((-(order as isize) + 1)..order as isize)
        .step_by(2)
        .map(|i| {
            let mi = F::from(i).unwrap();
            let num = unsafe { F::from(PI).unwrap_unchecked() * mi };
            let denom = unsafe { F::from(2. * order as f64).unwrap_unchecked() };
            let c = Complex::new(F::zero(), num / denom);
            -c.exp()
        })
        .collect::<Vec<_>>();

    ZpkFormatFilter::new(Vec::new(), p, F::one())
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[cfg(feature = "unstable")]
    #[test]
    fn matches_scipy_buttap() {
        let p: [Complex<f64>; 4] = [
            Complex::new(-0.38268343, 0.92387953),
            Complex::new(-0.92387953, 0.38268343),
            Complex::new(-0.92387953, -0.38268343),
            Complex::new(-0.38268343, -0.92387953),
        ];
        let zpk = buttap_st::<f64, 4>();
        for i in 0..4 {
            assert_relative_eq!(p[i].re, zpk.p[i].re, max_relative = 1e-7);
            assert_relative_eq!(p[i].im, zpk.p[i].im, max_relative = 1e-7);
        }
    }

    #[cfg(feature = "unstable")]
    #[test]
    fn matches_scipy_iirfilter_butter_zpk() {
        let expected_zpk: ZpkFormatFilter<f64, 8> = ZpkFormatFilter::new(
            Vec::from_slice(&[
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
            ])
            .unwrap(),
            Vec::from_slice(&[
                Complex::new(0.98924866, -0.03710237),
                Complex::new(0.96189799, -0.03364097),
                Complex::new(0.96189799, 0.03364097),
                Complex::new(0.98924866, 0.03710237),
                Complex::new(0.93873849, 0.16792939),
                Complex::new(0.89956011, 0.08396115),
                Complex::new(0.89956011, -0.08396115),
                Complex::new(0.93873849, -0.16792939),
            ])
            .unwrap(),
            2.6775767382597835e-5,
        );
        let filter = iirfilter_st::<f64, 4, _>(
            [10., 50.],
            None,
            None,
            Some(FilterBandType::Bandpass),
            Some(FilterType::Butterworth),
            Some(false),
            Some(FilterOutputType::Zpk),
            Some(1666.),
        );

        match filter {
            DigitalFilter::Zpk(zpk) => {
                assert_eq!(zpk.z.len(), expected_zpk.z.len());
                for (a, e) in zpk.z.iter().zip(expected_zpk.z.iter()) {
                    assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
                    assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
                }

                assert_eq!(zpk.p.len(), expected_zpk.p.len());
                for (a, e) in zpk.p.iter().zip(expected_zpk.p.iter()) {
                    assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
                    assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
                }

                assert_relative_eq!(zpk.k, expected_zpk.k, max_relative = 1e-8);
            }
            _ => panic!(),
        }
    }

    #[cfg(feature = "unstable")]
    #[test]
    fn matches_scipy_iirfilter_butter_sos() {
        let filter = iirfilter_st::<f64, 4, _>(
            [10., 50.],
            None,
            None,
            Some(FilterBandType::Bandpass),
            Some(FilterType::Butterworth),
            Some(false),
            Some(FilterOutputType::Sos),
            Some(1666.),
        );

        match filter {
            DigitalFilter::Sos(sos) => {
                // println!("{:?}", sos);

                let expected_sos = [
                    Sos::new(
                        [2.67757674e-05, 5.35515348e-05, 2.67757674e-05],
                        [1.00000000e+00, -1.79912022e+00, 8.16257861e-01],
                    ),
                    Sos::new(
                        [1.00000000e+00, 2.00000000e+00, 1.00000000e+00],
                        [1.00000000e+00, -1.87747699e+00, 9.09430241e-01],
                    ),
                    Sos::new(
                        [1.00000000e+00, -2.00000000e+00, 1.00000000e+00],
                        [1.00000000e+00, -1.92379599e+00, 9.26379467e-01],
                    ),
                    Sos::new(
                        [1.00000000e+00, -2.00000000e+00, 1.00000000e+00],
                        [1.00000000e+00, -1.97849731e+00, 9.79989489e-01],
                    ),
                ];

                assert_eq!(expected_sos.len(), sos.sos.len());
                for i in 0..sos.sos.len() {
                    let actual = sos.sos[i];
                    let expected = expected_sos[i];
                    assert_relative_eq!(actual.b[0], expected.b[0], max_relative = 1e-7);
                    assert_relative_eq!(actual.b[1], expected.b[1], max_relative = 1e-7);
                    assert_relative_eq!(actual.b[2], expected.b[2], max_relative = 1e-7);
                    assert_relative_eq!(actual.a[0], expected.a[0], max_relative = 1e-7);
                    assert_relative_eq!(actual.a[1], expected.a[1], max_relative = 1e-7);
                    assert_relative_eq!(actual.a[2], expected.a[2], max_relative = 1e-7);
                }
            }
            _ => panic!(),
        }
    }

    #[cfg(feature = "unstable")]
    #[test]
    fn matches_scipy_iirfilter_butter_ba() {
        let filter = iirfilter_st::<f64, 4, _>(
            [10., 50.],
            None,
            None,
            Some(FilterBandType::Bandpass),
            Some(FilterType::Butterworth),
            Some(false),
            Some(FilterOutputType::Ba),
            Some(1666.),
        );

        match filter {
            DigitalFilter::Ba(ba) => {
                let expected_b = [
                    2.67757674e-05,
                    0.00000000e+00,
                    -1.07103070e-04,
                    0.00000000e+00,
                    1.60654604e-04,
                    0.00000000e+00,
                    -1.07103070e-04,
                    0.00000000e+00,
                    2.67757674e-05,
                ];
                let expected_a = [
                    1.,
                    -7.57889051,
                    25.1632497,
                    -47.80506049,
                    56.83958432,
                    -43.31144279,
                    20.65538731,
                    -5.63674562,
                    0.67391808,
                ];

                assert_eq!(expected_b.len(), ba.b.len());
                assert_eq!(expected_a.len(), ba.a.len());
                for i in 0..expected_b.len() {
                    assert_relative_eq!(ba.b[0], expected_b[0], max_relative = 1e-7);
                    assert_relative_eq!(ba.b[1], expected_b[1], max_relative = 1e-7);
                    assert_relative_eq!(ba.b[2], expected_b[2], max_relative = 1e-7);
                    assert_relative_eq!(ba.b[3], expected_b[3], max_relative = 1e-7);
                    assert_relative_eq!(ba.b[4], expected_b[4], max_relative = 1e-7);

                    assert_relative_eq!(ba.a[0], expected_a[0], max_relative = 1e-7);
                    assert_relative_eq!(ba.a[1], expected_a[1], max_relative = 1e-7);
                    assert_relative_eq!(ba.a[2], expected_a[2], max_relative = 1e-7);
                    assert_relative_eq!(ba.a[3], expected_a[3], max_relative = 1e-7);
                    assert_relative_eq!(ba.a[4], expected_a[4], max_relative = 1e-7);
                }
            }
            _ => panic!(),
        }
    }

    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_buttap() {
        let p: [Complex<f64>; 4] = [
            Complex::new(-0.38268343, 0.92387953),
            Complex::new(-0.92387953, 0.38268343),
            Complex::new(-0.92387953, -0.38268343),
            Complex::new(-0.38268343, -0.92387953),
        ];
        let zpk = buttap_dyn::<f64>(4);
        for i in 0..4 {
            assert_relative_eq!(p[i].re, zpk.p[i].re, max_relative = 1e-7);
            assert_relative_eq!(p[i].im, zpk.p[i].im, max_relative = 1e-7);
        }
    }

    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_iirfilter_butter_zpk() {
        let expected_zpk: ZpkFormatFilter<f64> = ZpkFormatFilter::new(
            vec![
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
            ],
            vec![
                Complex::new(0.98924866, -0.03710237),
                Complex::new(0.96189799, -0.03364097),
                Complex::new(0.96189799, 0.03364097),
                Complex::new(0.98924866, 0.03710237),
                Complex::new(0.93873849, 0.16792939),
                Complex::new(0.89956011, 0.08396115),
                Complex::new(0.89956011, -0.08396115),
                Complex::new(0.93873849, -0.16792939),
            ],
            2.6775767382597835e-5,
        );
        let filter = iirfilter_dyn::<f64>(
            4,
            vec![10., 50.],
            None,
            None,
            Some(FilterBandType::Bandpass),
            Some(FilterType::Butterworth),
            Some(false),
            Some(FilterOutputType::Zpk),
            Some(1666.),
        );

        match filter {
            DigitalFilter::Zpk(zpk) => {
                assert_eq!(zpk.z.len(), expected_zpk.z.len());
                for (a, e) in zpk.z.iter().zip(expected_zpk.z.iter()) {
                    assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
                    assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
                }

                assert_eq!(zpk.p.len(), expected_zpk.p.len());
                for (a, e) in zpk.p.iter().zip(expected_zpk.p.iter()) {
                    assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
                    assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
                }

                assert_relative_eq!(zpk.k, expected_zpk.k, max_relative = 1e-8);
            }
            _ => panic!(),
        }
    }

    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_iirfilter_butter_sos() {
        let filter = iirfilter_dyn::<f64>(
            4,
            vec![10., 50.],
            None,
            None,
            Some(FilterBandType::Bandpass),
            Some(FilterType::Butterworth),
            Some(false),
            Some(FilterOutputType::Sos),
            Some(1666.),
        );

        match filter {
            DigitalFilter::Sos(sos) => {
                // println!("{:?}", sos);

                let expected_sos = [
                    Sos::new(
                        [2.67757674e-05, 5.35515348e-05, 2.67757674e-05],
                        [1.00000000e+00, -1.79912022e+00, 8.16257861e-01],
                    ),
                    Sos::new(
                        [1.00000000e+00, 2.00000000e+00, 1.00000000e+00],
                        [1.00000000e+00, -1.87747699e+00, 9.09430241e-01],
                    ),
                    Sos::new(
                        [1.00000000e+00, -2.00000000e+00, 1.00000000e+00],
                        [1.00000000e+00, -1.92379599e+00, 9.26379467e-01],
                    ),
                    Sos::new(
                        [1.00000000e+00, -2.00000000e+00, 1.00000000e+00],
                        [1.00000000e+00, -1.97849731e+00, 9.79989489e-01],
                    ),
                ];

                assert_eq!(expected_sos.len(), sos.sos.len());
                for i in 0..sos.sos.len() {
                    let actual = sos.sos[i];
                    let expected = expected_sos[i];
                    assert_relative_eq!(actual.b[0], expected.b[0], max_relative = 1e-7);
                    assert_relative_eq!(actual.b[1], expected.b[1], max_relative = 1e-7);
                    assert_relative_eq!(actual.b[2], expected.b[2], max_relative = 1e-7);
                    assert_relative_eq!(actual.a[0], expected.a[0], max_relative = 1e-7);
                    assert_relative_eq!(actual.a[1], expected.a[1], max_relative = 1e-7);
                    assert_relative_eq!(actual.a[2], expected.a[2], max_relative = 1e-7);
                }
            }
            _ => panic!(),
        }
    }

    #[cfg(feature = "use_std")]
    #[test]
    fn matches_scipy_iirfilter_butter_ba() {
        let filter = iirfilter_dyn::<f64>(
            4,
            vec![10., 50.],
            None,
            None,
            Some(FilterBandType::Bandpass),
            Some(FilterType::Butterworth),
            Some(false),
            Some(FilterOutputType::Ba),
            Some(1666.),
        );

        match filter {
            DigitalFilter::Ba(ba) => {
                let expected_b = [
                    2.67757674e-05,
                    0.00000000e+00,
                    -1.07103070e-04,
                    0.00000000e+00,
                    1.60654604e-04,
                    0.00000000e+00,
                    -1.07103070e-04,
                    0.00000000e+00,
                    2.67757674e-05,
                ];
                let expected_a = [
                    1.,
                    -7.57889051,
                    25.1632497,
                    -47.80506049,
                    56.83958432,
                    -43.31144279,
                    20.65538731,
                    -5.63674562,
                    0.67391808,
                ];

                assert_eq!(expected_b.len(), ba.b.len());
                assert_eq!(expected_a.len(), ba.a.len());
                assert_relative_eq!(ba.b[0], expected_b[0], max_relative = 1e-7);
                assert_relative_eq!(ba.b[1], expected_b[1], max_relative = 1e-7);
                assert_relative_eq!(ba.b[2], expected_b[2], max_relative = 1e-7);
                assert_relative_eq!(ba.b[3], expected_b[3], max_relative = 1e-7);
                assert_relative_eq!(ba.b[4], expected_b[4], max_relative = 1e-7);

                assert_relative_eq!(ba.a[0], expected_a[0], max_relative = 1e-7);
                assert_relative_eq!(ba.a[1], expected_a[1], max_relative = 1e-7);
                assert_relative_eq!(ba.a[2], expected_a[2], max_relative = 1e-7);
                assert_relative_eq!(ba.a[3], expected_a[3], max_relative = 1e-7);
                assert_relative_eq!(ba.a[4], expected_a[4], max_relative = 1e-7);
            }
            _ => panic!(),
        }
    }
}
