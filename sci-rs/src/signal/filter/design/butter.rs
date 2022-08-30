use nalgebra::RealField;
use num_traits::Float;

// use super::{iirfilter_st, FilterBandType, FilterOutputType, FilterType, Sos};

// ///
// /// Butterworth digital and analog filter design.
// ///
// /// Design an Nth-order digital or analog Butterworth filter and return
// /// the filter coefficients.
// ///
// /// Parameters
// /// ----------
// /// N : int
// ///     The order of the filter. For 'bandpass' and 'bandstop' filters,
// ///     the resulting order of the final second-order sections ('sos')
// ///     matrix is ``2*N``, with `N` the number of biquad sections
// ///     of the desired system.
// /// Wn : array_like
// ///     The critical frequency or frequencies. For lowpass and highpass
// ///     filters, Wn is a scalar; for bandpass and bandstop filters,
// ///     Wn is a length-2 sequence.
// ///
// ///     For a Butterworth filter, this is the point at which the gain
// ///     drops to 1/sqrt(2) that of the passband (the "-3 dB point").
// ///
// ///     For digital filters, if `fs` is not specified, `Wn` units are
// ///     normalized from 0 to 1, where 1 is the Nyquist frequency (`Wn` is
// ///     thus in half cycles / sample and defined as 2*critical frequencies
// ///     / `fs`). If `fs` is specified, `Wn` is in the same units as `fs`.
// ///
// ///     For analog filters, `Wn` is an angular frequency (e.g. rad/s).
// /// btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
// ///     The type of filter.  Default is 'lowpass'.
// /// analog : bool, optional
// ///     When True, return an analog filter, otherwise a digital filter is
// ///     returned.
// /// output : {'ba', 'zpk', 'sos'}, optional
// ///     Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
// ///     second-order sections ('sos'). Default is 'ba' for backwards
// ///     compatibility, but 'sos' should be used for general-purpose filtering.
// /// fs : float, optional
// ///     The sampling frequency of the digital system.
// ///
// pub fn butter_st<F, const N: usize, const M: usize>(
//     n: usize,
//     wn: [F; M],
//     btype: Option<FilterBandType>,
//     analog: Option<bool>,
//     output: Option<FilterOutputType>,
//     fs: Option<F>,
// ) -> FilterOutput<F, N, { N * 2 }>
// where
//     F: RealField + Float,
// {
//     let btype = btype.unwrap_or(FilterBandType::Lowpass);
//     let analog = analog.unwrap_or(false);
//     let output = output.unwrap_or(FilterOutputType::Ba);
//     iirfilter_st::<F, N, M>(
//         wn,
//         None,
//         None,
//         Some(btype),
//         Some(FilterType::Butterworth),
//         Some(analog),
//         Some(output),
//         fs,
//     )
// }
