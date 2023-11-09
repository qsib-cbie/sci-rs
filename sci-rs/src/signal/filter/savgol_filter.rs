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
pub fn savgol_filter_dyn<YI, F>(y: YI, window_length: usize, polyorder: usize) -> Vec<F>
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

    let fir = savgol_coeffs_dyn::<f64>(window_length, polyorder)
        .into_iter()
        .map(|f| F::from_f64(f).unwrap())
        .collect::<Vec<_>>();

    // Pad with nearest edge value
    let size_hint = y.size_hint();
    let mut data = Vec::with_capacity(size_hint.1.unwrap_or(size_hint.0));
    let mut y = y;
    let Some(nearest) = y.next() else {
        return vec![];
    };
    let nearest = *nearest.borrow();
    data.extend((0..window_length / 2).map(|_| nearest));
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
pub fn savgol_coeffs_dyn<F>(window_length: usize, polyorder: usize) -> Vec<F>
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
            .map(|i| (F::from_usize(i).unwrap() + f - half_window))
            .collect::<Vec<_>>()
    } else {
        (0..window_length)
            .map(|i| F::from_usize(i).unwrap() - half_window)
            .collect::<Vec<_>>()
    };

    // Columns are 2m+1 integer positions centered on 0
    // Rows are powers of positions from 0 to polyorder
    // Setting up a Vandermonde matrix for solving A * coeffs = y
    #[allow(non_snake_case)]
    let A = na::DMatrix::<F>::from_fn(polyorder + 1, window_length, |i, j| pos[j].powi(i as i32));
    let y = na::DVector::<F>::from_fn(
        polyorder + 1,
        |i, _| {
            if i == 0 {
                F::one()
            } else {
                F::zero()
            }
        },
    );

    // Solve the system for the Savitsky-Golay FIR coefficients
    let solve = lstsq::lstsq(&A, &y, F::from_f32(1e-9).unwrap()).unwrap();
    solve.solution.data.into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[test]
    pub fn can_filter() {
        let v = savgol_filter_dyn((0..100).map(|i| i as f32), 11, 2);
        println!("v = {:?}", v);

        let v = savgol_filter_dyn((0..0).map(|i| i as f32), 11, 2);
        println!("v = {:?}", v);
    }

    #[test]
    pub fn can_coeffs() {
        let actual = savgol_coeffs_dyn::<f32>(5, 2);
        let expected = [-0.08571429, 0.34285714, 0.48571429, 0.34285714, -0.08571429];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-5);
        }

        let actual = savgol_coeffs_dyn::<f64>(5, 2);
        let expected = [-0.08571429, 0.34285714, 0.48571429, 0.34285714, -0.08571429];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-7);
        }

        let actual = savgol_coeffs_dyn::<f64>(4, 2);
        let expected = [-0.0625, 0.5625, 0.5625, -0.0625];
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a, e, max_relative = 1e-7);
        }

        let actual = savgol_coeffs_dyn::<f64>(21, 8);
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
    }
}
