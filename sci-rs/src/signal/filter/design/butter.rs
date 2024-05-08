use core::iter::Sum;

use nalgebra::RealField;
use num_traits::Float;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

#[cfg(feature = "alloc")]
use super::{iirfilter_dyn, DigitalFilter, FilterBandType, FilterOutputType, FilterType, Sos};
///
/// Butterworth digital and analog filter design.
///
/// Design an Nth-order digital or analog Butterworth filter and return
/// the filter coefficients.
///
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html>
///
#[cfg(feature = "alloc")]
pub fn butter_dyn<F>(
    order: usize,
    wn: Vec<F>,
    btype: Option<FilterBandType>,
    analog: Option<bool>,
    output: Option<FilterOutputType>,
    fs: Option<F>,
) -> DigitalFilter<F>
where
    F: RealField + Float + Sum,
{
    let btype = btype.unwrap_or(FilterBandType::Lowpass);
    let analog = analog.unwrap_or(false);
    let output = output.unwrap_or(FilterOutputType::Ba);
    iirfilter_dyn::<F>(
        order,
        wn,
        None,
        None,
        Some(btype),
        Some(FilterType::Butterworth),
        Some(analog),
        Some(output),
        fs,
    )
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Complex;

    use crate::signal::filter::design::ZpkFormatFilter;

    use super::*;

    #[cfg(feature = "alloc")]
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
        let filter = butter_dyn::<f64>(
            4,
            vec![10., 50.],
            Some(FilterBandType::Bandpass),
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

    #[cfg(feature = "alloc")]
    #[test]
    fn matches_scipy_iirfilter_butter_sos() {
        let filter = butter_dyn::<f64>(
            4,
            vec![10., 50.],
            Some(FilterBandType::Bandpass),
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

    #[cfg(feature = "alloc")]
    #[test]
    fn matches_scipy_iirfilter_butter_sos_highpass() {
        let filter = butter_dyn::<f64>(
            5,
            vec![2.],
            Some(FilterBandType::Highpass),
            Some(false),
            Some(FilterOutputType::Sos),
            Some(600.),
        );
        // >>> butter(5,[2],btype='hp',output='sos',fs=600)
        // array([[ 0.966679  , -0.966679  ,  0.        ,  1.        , -0.97927235, 0.        ],
        //        [ 1.        , -2.        ,  1.        ,  1.        , -1.96624768, 0.966679  ],
        //        [ 1.        , -2.        ,  1.        ,  1.        , -1.98670428, 0.9871401 ]])
        match filter {
            DigitalFilter::Sos(sos) => {
                // println!("{:?}", sos);

                let expected_sos = [
                    Sos::new([0.966679, -0.966679, 0.], [1., -0.97927235, 0.]),
                    Sos::new([1., -2., 1.], [1., -1.96624768, 0.966679]),
                    Sos::new([1., -2., 1.], [1., -1.98670428, 0.9871401]),
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

    #[cfg(feature = "alloc")]
    #[test]
    fn matches_scipy_iirfilter_butter_ba() {
        let filter = butter_dyn::<f64>(
            4,
            vec![10., 50.],
            Some(FilterBandType::Bandpass),
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
