use core::{f64::consts::PI, ops::Mul};

use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{
    bilinear_zpk, lp2bp_zpk, zpk2sos, FilterBandType, FilterOutput, FilterOutputType, FilterType,
    Sos, Zpk,
};

// filter_dict = {'butter': [buttap, buttord],
//                'butterworth': [buttap, buttord],
//                'cauer': [ellipap, ellipord],
//                'elliptic': [ellipap, ellipord],
//                'ellip': [ellipap, ellipord],
//                'bessel': [besselap],
//                'bessel_phase': [besselap],
//                'bessel_delay': [besselap],
//                'bessel_mag': [besselap],
//                'cheby1': [cheb1ap, cheb1ord],
//                'chebyshev1': [cheb1ap, cheb1ord],
//                'chebyshevi': [cheb1ap, cheb1ord],
//                'cheby2': [cheb2ap, cheb2ord],
//                'chebyshev2': [cheb2ap, cheb2ord],
//                'chebyshevii': [cheb2ap, cheb2ord],
//                }

///
///
/// IIR digital and analog filter design given order and critical points.
///
/// Design an Nth-order digital or analog filter and return the filter
/// coefficients.
///
/// Parameters
/// ----------
/// N : int
///     The order of the filter.
/// Wn : array_like
///     A scalar or length-2 sequence giving the critical frequencies.
///
///     For digital filters, `Wn` are in the same units as `fs`. By default,
///     `fs` is 2 half-cycles/sample, so these are normalized from 0 to 1,
///     where 1 is the Nyquist frequency. (`Wn` is thus in
///     half-cycles / sample.)
///
///     For analog filters, `Wn` is an angular frequency (e.g., rad/s).
///
///     When Wn is a length-2 sequence, ``Wn[0]`` must be less than ``Wn[1]``.
/// rp : float, optional
///     For Chebyshev and elliptic filters, provides the maximum ripple
///     in the passband. (dB)
/// rs : float, optional
///     For Chebyshev and elliptic filters, provides the minimum attenuation
///     in the stop band. (dB)
/// btype : {'bandpass', 'lowpass', 'highpass', 'bandstop'}, optional
///     The type of filter.  Default is 'bandpass'.
/// analog : bool, optional
///     When True, return an analog filter, otherwise a digital filter is
///     returned.
/// ftype : str, optional
///     The type of IIR filter to design:
///
///         - Butterworth   : 'butter'
///         - Chebyshev I   : 'cheby1'
///         - Chebyshev II  : 'cheby2'
///         - Cauer/elliptic: 'ellip'
///         - Bessel/Thomson: 'bessel'
///
/// output : {'ba', 'zpk', 'sos'}, optional
///     Filter form of the output:
///
///         - second-order sections (recommended): 'sos'
///         - numerator/denominator (default)    : 'ba'
///         - pole-zero                          : 'zpk'
///
///     In general the second-order sections ('sos') form  is
///     recommended because inferring the coefficients for the
///     numerator/denominator form ('ba') suffers from numerical
///     instabilities. For reasons of backward compatibility the default
///     form is the numerator/denominator form ('ba'), where the 'b'
///     and the 'a' in 'ba' refer to the commonly used names of the
///     coefficients used.
///
///     Note: Using the second-order sections form ('sos') is sometimes
///     associated with additional computational costs: for
///     data-intense use cases it is therefore recommended to also
///     investigate the numerator/denominator form ('ba').
///
/// fs : float, optional
///     The sampling frequency of the digital system.
///
///     .. versionadded:: 1.2.0
///
/// Returns
/// -------
/// b, a : ndarray, ndarray
///     Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
///     Only returned if ``output='ba'``.
/// z, p, k : ndarray, ndarray, float
///     Zeros, poles, and system gain of the IIR filter transfer
///     function.  Only returned if ``output='zpk'``.
/// sos : ndarray
///     Second-order sections representation of the IIR filter.
///     Only returned if ``output=='sos'``.
///
/// See Also
/// --------
/// butter : Filter design using order and critical points
/// cheby1, cheby2, ellip, bessel
/// buttord : Find order and critical points from passband and stopband spec
/// cheb1ord, cheb2ord, ellipord
/// iirdesign : General filter design using passband and stopband spec
///
pub fn iirfilter_st<F, const N: usize, const M: usize>(
    wn: [F; M],
    rp: Option<F>,
    rs: Option<F>,
    btype: Option<FilterBandType>,
    ftype: Option<FilterType>,
    analog: Option<bool>,
    output: Option<FilterOutputType>,
    fs: Option<F>,
) -> FilterOutput<F, N>
where
    F: RealField + Float,
    [Sos<F>; N / 2 - 1]: Sized,
{
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

    if M > 1 && wn[0] < wn[1] {
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
    let zpk = match ftype {
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
            lp2bp_zpk(zpk, Some(wo), Some(bw))
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
        FilterOutputType::Zpk => FilterOutput::Zpk(zpk),
        FilterOutputType::Ba => {
            // zpk2tf(z, p, k),
            todo!()
        }
        FilterOutputType::Sos => FilterOutput::Sos(zpk2sos(zpk, analog)),
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
fn buttap<F, const N: usize>() -> Zpk<F, N>
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

    Zpk::new(heapless::Vec::new(), p, F::one())
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
}
