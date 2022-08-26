use core::{f64::consts::PI, ops::Mul};

use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use crate::signal::filter::design::ZpkFormatFilter;

use super::{
    bilinear_zpk, lp2bp_zpk, zpk2sos, DigitalFilter, FilterBandType, FilterOutputType, FilterType,
    Sos,
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
pub fn iirfilter_st<F, const N: usize, const N2: usize, const M: usize>(
    wn: [F; M],
    rp: Option<F>,
    rs: Option<F>,
    btype: Option<FilterBandType>,
    ftype: Option<FilterType>,
    analog: Option<bool>,
    output: Option<FilterOutputType>,
    fs: Option<F>,
) -> DigitalFilter<F, N, N2>
where
    F: RealField + Float,
{
    assert!(N * 2 == N2);

    let analog = analog.unwrap_or(false);
    let mut wn = wn;

    if let Some(fs) = fs {
        if analog {
            panic!("fs cannot be specified for an analog filter");
        }

        for i in 0..M {
            wn[i] = F::from(2.).unwrap() * wn[i] / fs;
        }
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
        FilterType::Butterworth => buttap::<F, N>(),
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
        (fs.unwrap_or(F::one()), wn)
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
            lp2bp_zpk::<F, N, N2>(zpk, Some(wo), Some(bw))
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
        FilterOutputType::Ba => {
            // zpk2tf(z, p, k),
            todo!()
        }
        FilterOutputType::Sos => DigitalFilter::Sos(zpk2sos(zpk, None, Some(analog))),
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
fn buttap<F, const N: usize>() -> ZpkFormatFilter<F, N>
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

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn matches_scipy_buttap() {
        let p: [Complex<f64>; 4] = [
            Complex::new(-0.38268343, 0.92387953),
            Complex::new(-0.92387953, 0.38268343),
            Complex::new(-0.92387953, -0.38268343),
            Complex::new(-0.38268343, -0.92387953),
        ];
        let zpk = buttap::<f64, 4>();
        for i in 0..4 {
            assert_relative_eq!(p[i].re, zpk.p[i].re, max_relative = 1e-7);
            assert_relative_eq!(p[i].im, zpk.p[i].im, max_relative = 1e-7);
        }
    }

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
        let filter = iirfilter_st::<f64, 4, 8, _>(
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

    // #[ignore]
    // #[test]
    // fn matches_scipy_iirfilter_butter_sos() {
    //     let filter = iirfilter_st::<f64, 4, 8, _>(
    //         [10., 50.],
    //         None,
    //         None,
    //         Some(FilterBandType::Bandpass),
    //         Some(FilterType::Butterworth),
    //         Some(false),
    //         Some(FilterOutputType::Sos),
    //         Some(1666.),
    //     );

    //     match filter {
    //         DigitalFilter::Sos(sos) => {}
    //         _ => panic!(),
    //     }
    // }
}
