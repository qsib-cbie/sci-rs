use core::iter::Sum;
use core::{
    borrow::Borrow,
    mem::{transmute, MaybeUninit},
};
use nalgebra as na;
use nalgebra::RealField;
use num_traits::Float;

#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

///
/// Savitzky-Golay filtering
///
/// <https://doi.org/10.1021/ac60214a047>
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter>
///
/// Design coefficients for a Savitzky-Golay filter and convolve with data using nearest edge padding.
///
pub fn savgol_filter_dyn<YI, F>(
    y: YI,
    window_length: usize,
    polyorder: usize,
    deriv: Option<usize>,
    delta: Option<F>,
) -> Vec<F>
where
    F: RealField + Copy + Sum,
    YI: Iterator,
    YI::Item: Borrow<F>,
{
    if window_length % 2 == 0 {
        panic!("window_length must be odd")
    }

    if window_length < polyorder + 2 {
        panic!("window_length is too small for the polynomials order")
    }

    

    let mut fir = savgol_coeffs_dyn::<F>(window_length, polyorder, deriv, delta)
        .into_iter()
        .collect::<Vec<_>>();

    fir.reverse();

    // Pad with nearest edge value
    let size_hint = y.size_hint();
    let mut data = Vec::with_capacity(size_hint.1.unwrap_or(size_hint.0));
    let mut y = y;
    let Some(nearest) = y.next() else {
        return vec![];
    };
    let nearest = *nearest.borrow();
    data.extend((0..(window_length / 2 + 1)).map(|_| nearest));
    data.extend(y.map(|yi| *yi.borrow()));
    let nearest = *data.last().unwrap();
    data.extend((0..window_length / 2).map(|_| nearest));

    // Convolve the data with the FIR coefficients
    let rslt = data
        .windows(window_length)
        .map(|w| w.iter().zip(fir.iter()).map(|(a, b)| *a * *b).sum::<F>())
        .collect();

    rslt
}

///
/// Design 1-D Savitzky-Golay filter coefficients
///
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_coeffs.html#scipy.signal.savgol_coeffs>
///
/// This function is sensitive to f64 and f32 primitives due to use of least squares.
/// The coefficients may go to zero for higher order polynomials and larger window lengths.
///
pub fn savgol_coeffs_dyn<F>(
    window_length: usize,
    polyorder: usize,
    deriv: Option<usize>,
    delta: Option<F>,
) -> Vec<F>
where
    F: RealField + Copy,
{
    if polyorder >= window_length {
        panic!("polyorder must be less than window_length")
    }

    let half_window = F::from_usize(window_length / 2).unwrap();
    let rem = window_length % 2;

    let pos = if rem == 0 {
        let f = F::from_f32(0.5).unwrap();
        (0..window_length)
            .map(|i| (half_window - F::from_usize(i).unwrap() - f))
            .collect::<Vec<_>>()
    } else {
        (0..window_length)
            .map(|i| half_window - F::from_usize(i).unwrap())
            .collect::<Vec<_>>()
    };

    //handle the case of default args
    let der = deriv.unwrap_or(0);
    let del = delta.unwrap_or(F::one());

    if der > polyorder {
        let mut ret = vec![F::zero(); window_length];
        return ret;
    }

    // Columns are 2m+1 integer positions centered on 0
    // Rows are powers of positions from 0 to polyorder
    // Setting up a Vandermonde matrix for solving A * coeffs = y
    #[allow(non_snake_case)]
    let A = na::DMatrix::<F>::from_fn(polyorder + 1, window_length, |i, j| pos[j].powi(i as i32));
    let mut y =
        na::DVector::<F>::from_fn(
            polyorder + 1,
            |i, _| {
                if i == der {
                    F::one()
                } else {
                    F::zero()
                }
            },
        );

    y[der] = (F::from_usize(factorial(der)).unwrap()) / del.powi(der as i32);

    // Solve the system for the Savitsky-Golay FIR coefficients
    let solve = lstsq::lstsq(&A, &y, F::from_f32(1e-9).unwrap()).unwrap();
    solve.solution.data.into()
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[cfg(feature = "std")]
    #[test]
    pub fn can_filter() {
        let v = savgol_filter_dyn((0..100).map(|i| i as f32), 11, 2, None, None);
        println!("v = {:?}", v);

        let v = savgol_filter_dyn((0..0).map(|i| i as f32), 11, 2, None, None);
        println!("v = {:?}", v);

        let actual_coeff = savgol_coeffs_dyn::<f64>(5, 2, Some(0), Some(1.0));
        println!("coeffs = {:?}", actual_coeff);
        let input = [2.0, 2.0, 5.0, 2.0, 1.0, 0.0, 1.0, 4.0, 9.0];
        let actual = savgol_filter_dyn(input.iter(), 5, 2, None, None);
        println!("actual = {:?}", actual);
        let expected = [
            1.74285714, 3.02857143, 3.54285714, 2.85714286, 0.65714286, 0.17142857, 1.0, 4.6,
            7.97142857,
        ];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-5);
        }

        let actual_coeff = savgol_coeffs_dyn::<f64>(5, 2, Some(1), None);
        println!("coeffs = {:?}", actual_coeff);
        let input = [2.0, 2.0, 5.0, 2.0, 1.0, 0.0, 1.0, 4.0, 9.0];
        let actual = savgol_filter_dyn(input.iter(), 5, 2, Some(1), None);
        println!("actual = {:?}", actual);
        let expected = [0.6, 0.3, -0.2, -0.8, -1.0, 0.4, 2.0, 2.6, 2.1];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-5);
        }

        let input = (0..100).map(|i| (3 * i - 2) as f64); //y = 3x - 2
        let actual = savgol_filter_dyn(input, 51, 5, None, None);
        let expected = [
            2.45650177,
            4.06008889,
            5.86936071,
            7.87987343,
            10.08430349,
            12.47257185,
            15.03201779,
            17.74762253,
            20.6022825,
            23.57713222,
            26.65191695,
            29.805415,
            33.01590968,
            36.26171099,
            39.52172696,
            42.77608466,
            46.00680095,
            49.19850285,
            52.33919758,
            55.4210924,
            58.44146398,
            61.40357756,
            64.31765571,
            67.20189688,
            70.08354354,
            73.,
            76.,
            79.,
            82.,
            85.,
            88.,
            91.,
            94.,
            97.,
            100.,
            103.,
            106.,
            109.,
            112.,
            115.,
            118.,
            121.,
            124.,
            127.,
            130.,
            133.,
            136.,
            139.,
            142.,
            145.,
            148.,
            151.,
            154.,
            157.,
            160.,
            163.,
            166.,
            169.,
            172.,
            175.,
            178.,
            181.,
            184.,
            187.,
            190.,
            193.,
            196.,
            199.,
            202.,
            205.,
            208.,
            211.,
            214.,
            217.,
            220.,
            222.91645647,
            225.79810312,
            228.6823443,
            231.59642245,
            234.55853602,
            237.5789076,
            240.66080242,
            243.80149716,
            246.99319905,
            250.22391534,
            253.47827305,
            256.73828901,
            259.98409032,
            263.194585,
            266.34808305,
            269.42286778,
            272.3977175,
            275.25237747,
            277.96798222,
            280.52742816,
            282.91569651,
            285.12012658,
            287.13063929,
            288.93991111,
            290.54349823,
        ];

        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-5);
        }
    }

    #[test]
    pub fn can_coeffs() {
        let actual = savgol_coeffs_dyn::<f32>(5, 2, None, None);
        let expected = [-0.08571429, 0.34285714, 0.48571429, 0.34285714, -0.08571429];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-5);
        }

        let actual = savgol_coeffs_dyn::<f64>(5, 2, Some(0), Some(1.0));
        let expected = [-0.08571429, 0.34285714, 0.48571429, 0.34285714, -0.08571429];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-7);
        }

        let actual = savgol_coeffs_dyn::<f64>(4, 2, None, None);
        let expected = [-0.0625, 0.5625, 0.5625, -0.0625];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-7);
        }

        let actual = savgol_coeffs_dyn::<f64>(51, 5, None, None);
        let expected = [
            0.02784785,
            0.01160327,
            -0.00086484,
            -0.00994566,
            -0.01601181,
            -0.01941934,
            -0.02050775,
            -0.01959997,
            -0.01700239,
            -0.0130048,
            -0.00788047,
            -0.00188609,
            0.00473822,
            0.01176888,
            0.01899888,
            0.02623777,
            0.03331167,
            0.04006325,
            0.04635174,
            0.05205294,
            0.05705919,
            0.06127943,
            0.06463912,
            0.0670803,
            0.06856157,
            0.06905808,
            0.06856157,
            0.0670803,
            0.06463912,
            0.06127943,
            0.05705919,
            0.05205294,
            0.04635174,
            0.04006325,
            0.03331167,
            0.02623777,
            0.01899888,
            0.01176888,
            0.00473822,
            -0.00188609,
            -0.00788047,
            -0.0130048,
            -0.01700239,
            -0.01959997,
            -0.02050775,
            -0.01941934,
            -0.01601181,
            -0.00994566,
            -0.00086484,
            0.01160327,
            0.02784785,
        ];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 5e-6);
        }

        let actual = savgol_coeffs_dyn::<f64>(21, 8, None, None);
        let expected = [
            0.0125937,
            -0.04897551,
            0.03811252,
            0.04592441,
            -0.01514104,
            -0.06782274,
            -0.05517056,
            0.03024958,
            0.15283999,
            0.25791748,
            0.29894434,
            0.25791748,
            0.15283999,
            0.03024958,
            -0.05517056,
            -0.06782274,
            -0.01514104,
            0.04592441,
            0.03811252,
            -0.04897551,
            0.0125937,
        ];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-6);
        }

        //deriv tests
        let actual = savgol_coeffs_dyn::<f64>(5, 2, Some(1), Some(1.0));
        let expected = [2.0e-1, 1.0e-1, 2.07548111e-16, -1.0e-1, -2.0e-1];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-7);
        }

        let actual = savgol_coeffs_dyn::<f64>(6, 3, Some(1), None);
        let expected = [
            -0.09093915,
            0.4130291,
            0.21560847,
            -0.21560847,
            -0.4130291,
            0.09093915,
        ];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-7);
        }

        let actual = savgol_coeffs_dyn::<f64>(6, 3, Some(2), Some(1.0));
        let expected = [
            0.17857143,
            -0.03571429,
            -0.14285714,
            -0.14285714,
            -0.03571429,
            0.17857143,
        ];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 5e-6);
        }
    }
}
