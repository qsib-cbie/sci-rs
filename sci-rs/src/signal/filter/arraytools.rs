//! Functions for acting on a axis of an array.
//!
//! Designed for ndarrays; with scipy's internal nomenclature.

use alloc::{vec, vec::Vec};
use ndarray::{
    ArrayBase, ArrayView, Axis, Data, Dim, Dimension, IntoDimension, Ix, RemoveAxis, SliceArg,
    SliceInfo, SliceInfoElem,
};
use sci_rs_core::{Error, Result};

/// Internal function for casting into [Axis] and appropriate usize from isize.
///
/// # Parameters
/// axis: The user-specificed axis which filter is to be applied on.
/// x: The input-data whose axis object that will be manipulated against.
///
/// # Notes
/// Const nature of this function means error has to be manually created.
#[inline]
pub(crate) const fn check_and_get_axis_st<'a, T, S, const N: usize>(
    axis: Option<isize>,
    x: &ArrayBase<S, Dim<[Ix; N]>>,
) -> core::result::Result<usize, ()>
where
    S: Data<Elem = T> + 'a,
{
    // Before we convert into the appropriate axis object, we have to check at runtime that the
    // axis value specified is within -N <= axis < N.
    match axis {
        None => (),
        Some(axis) if axis.is_negative() && axis.unsigned_abs() > N => {
            return Err(());
        }
        Some(axis) if axis.is_negative() && axis.unsigned_abs() <= N => (),
        Some(axis) if axis.unsigned_abs() >= N => {
            return Err(());
        }
        Some(axis) if axis.unsigned_abs() < N => (),
        Some(axis) => panic!("This should not be reached!"),
    }

    // We make a best effort to convert into appropriate axis object.
    let axis_inner: isize = match axis {
        Some(axis) => axis,
        None => -1,
    };
    if axis_inner >= 0 {
        Ok(axis_inner.unsigned_abs())
    } else {
        let axis_inner = N
            .checked_add_signed(axis_inner)
            .expect("Invalid add to `axis` option");
        Ok(axis_inner)
    }
}

/// Internal function for casting into [Axis] and appropriate usize from isize.
/// [check_and_get_axis_st] but without const, especially for IxDyn arrays.
///
/// # Parameters
/// axis: The user-specificed axis which filter is to be applied on.
/// x: The input-data whose axis object that will be manipulated against.
#[inline]
pub(crate) fn check_and_get_axis_dyn<'a, T, S, D>(
    axis: Option<isize>,
    x: &ArrayBase<S, D>,
) -> Result<usize>
where
    D: Dimension,
    S: Data<Elem = T> + 'a,
{
    let ndim = D::NDIM.unwrap_or(x.ndim());
    // Before we convert into the appropriate axis object, we have to check at runtime that the
    // axis value specified is within -N <= axis < N.
    if axis.is_some_and(|axis| {
        !(if axis < 0 {
            axis.unsigned_abs() <= ndim
        } else {
            axis.unsigned_abs() < ndim
        })
    }) {
        return Err(Error::InvalidArg {
            arg: "axis".into(),
            reason: "index out of range.".into(),
        });
    }

    // We make a best effort to convert into appropriate axis object.
    let axis_inner: isize = axis.unwrap_or(-1);
    if axis_inner >= 0 {
        Ok(axis_inner.unsigned_abs())
    } else {
        let axis_inner = ndim
            .checked_add_signed(axis_inner)
            .expect("Invalid add to `axis` option");
        Ok(axis_inner)
    }
}

/// Internal function for obtaining length of all axis as array from input from input.
///
/// This is almost the same as `a.shape()`, but is a array `[T; N]` instead of a slice `&[T]`.
///
/// # Parameters
/// `a`: Array whose shape is needed as a slice.
pub(crate) fn ndarray_shape_as_array_st<'a, S, T, const N: usize>(
    a: &ArrayBase<S, Dim<[Ix; N]>>,
) -> [Ix; N]
where
    [Ix; N]: IntoDimension<Dim = Dim<[Ix; N]>>,
    Dim<[Ix; N]>: RemoveAxis,
    S: Data<Elem = T> + 'a,
{
    a.shape().try_into().expect("Could not cast shape to array")
}

/// Takes a slice along `axis` from `a`.
///
/// # Parameters
/// * `a`: Array being sliced from.
/// * `start`: `Option<isize>`. None defaults to 0.
/// * `end`: `Option<isize>`.
/// * `step`: `Option<isize>`. None default to 1.
/// * `axis`: `Option<isize>`. None defaults to -1.
///
/// # Errors
/// - Axis is out of bounds.
///
/// # Panics
/// - Start/stop elements are out of bounds.
pub fn axis_slice<A, S, D>(
    a: &ArrayBase<S, D>,
    start: Option<isize>,
    end: Option<isize>,
    step: Option<isize>,
    axis: Option<isize>,
) -> Result<ArrayView<'_, A, D>>
where
    S: Data<Elem = A>,
    D: Dimension,
    SliceInfo<Vec<SliceInfoElem>, D, D>: SliceArg<D, OutDim = D>,
{
    let ndim = D::NDIM.unwrap_or(a.ndim());

    let axis = {
        if axis.is_some_and(|axis| {
            !(if axis < 0 {
                axis.unsigned_abs() <= ndim
            } else {
                axis.unsigned_abs() < ndim
            })
        }) {
            return Err(Error::InvalidArg {
                arg: "axis".into(),
                reason: "index out of range.".into(),
            });
        }

        // We make a best effort to convert into appropriate usize.
        let axis: isize = axis.unwrap_or(-1);
        if axis >= 0 {
            axis.unsigned_abs()
        } else {
            a.ndim()
                .checked_add_signed(axis)
                .expect("Invalid add to `axis` option")
        }
    };

    unsafe { axis_slice_unsafe(a, start, end, step, axis, ndim) }
}

/// Takes a slice along `axis` from `a`.
///
/// Assumes that the specified axis is within bounds.
///
/// # Parameters
/// * `a`: Array being sliced from.
/// * `start`: `Option<isize>`. None defaults to 0.
/// * `end`: `Option<isize>`.
/// * `step`: `Option<isize>`. None default to 1.
/// * `axis`: `usize`.
/// * `a_ndim`: Dimensionality of `a`. This strictly has to be `a.ndim()`.
///
/// # Panics
/// - Axis is out of bounds.
/// - Start/stop elements are out of bounds.
pub(crate) fn axis_slice_unsafe<A, S, D>(
    a: &ArrayBase<S, D>,
    start: Option<isize>,
    end: Option<isize>,
    step: Option<isize>,
    axis: usize,
    a_ndim: usize,
) -> Result<ArrayView<'_, A, D>>
where
    S: Data<Elem = A>,
    D: Dimension,
    SliceInfo<Vec<SliceInfoElem>, D, D>: SliceArg<D, OutDim = D>,
{
    debug_assert!(if a_ndim != 0 {
        axis < a_ndim // Eg: A 1D-array should only have axis = 0
    } else {
        axis <= a_ndim // Allow for axis = 0 when ndim = 0.
    });

    let axis_len = a.shape()[axis] as isize;
    let step = step.unwrap_or(1);

    let coerce = |idx: Option<isize>, def_pos: isize, def_neg: isize| -> isize {
        match idx {
            Some(i) if i.is_negative() => (axis_len + i),
            Some(i) => i.min(axis_len),
            None => {
                if !step.is_negative() {
                    def_pos
                } else {
                    def_neg
                }
            }
        }
    };
    let (start, end) = {
        let mut start = coerce(start, 0, axis_len - 1);
        let mut end = coerce(end, axis_len, -1);
        if step.is_negative() {
            (end + 1, Some(start + 1))
        } else {
            (start, Some(end)) // No + 1 breaking into axis_len
        }
    };

    let sl = SliceInfo::<_, D, D>::try_from({
        let mut tmp = vec![SliceInfoElem::from(..); a_ndim];
        tmp[axis] = SliceInfoElem::Slice { start, end, step };

        tmp
    })
    .unwrap();

    Ok(a.slice(&sl))
}

/// Reverse the 1-D slices (aka lanes) of `a` along axis `axis`.
///
/// Returns axis_slice(a, step=-1, axis=axis).
///
/// # Parameters
/// * `a`: Array being sliced from.
/// * `axis`: `Option<isize>`. None defaults to -1.
pub fn axis_reverse<A, S, D>(
    a: &ArrayBase<S, D>,
    axis: Option<isize>,
) -> Result<ArrayView<'_, A, D>>
where
    S: Data<Elem = A>,
    D: Dimension,
    SliceInfo<Vec<SliceInfoElem>, D, D>: SliceArg<D, OutDim = D>,
{
    let ndim = D::NDIM.unwrap_or(a.ndim());

    let axis = {
        if axis.is_some_and(|axis| {
            !(if axis < 0 {
                axis.unsigned_abs() <= ndim
            } else {
                axis.unsigned_abs() < ndim
            })
        }) {
            return Err(Error::InvalidArg {
                arg: "axis".into(),
                reason: "index out of range.".into(),
            });
        }

        // We make a best effort to convert into appropriate usize.
        let axis: isize = axis.unwrap_or(-1);
        if axis >= 0 {
            axis.unsigned_abs()
        } else {
            a.ndim()
                .checked_add_signed(axis)
                .expect("Invalid add to `axis` option")
        }
    };

    unsafe { axis_slice_unsafe(a, None, None, Some(-1), axis, ndim) }
}

/// Reverse the 1-D slices (aka lanes) of `a` along axis `axis`.
///
/// Returns axis_slice(a, step=-1, axis=axis).
///
/// # Parameters
/// * `a`: Array being sliced from.
/// * `axis`: `usize`.
/// * `a_ndim`: Dimensionality of `a`. This strictly has to be `a.ndim()`.
///
/// # Panics
/// If axis is out of bounds, and dimensions are wrong.
#[inline]
pub(crate) unsafe fn axis_reverse_unsafe<A, S, D>(
    a: &ArrayBase<S, D>,
    axis: usize,
    a_ndim: usize,
) -> ArrayView<'_, A, D>
where
    S: Data<Elem = A>,
    D: Dimension,
    SliceInfo<Vec<SliceInfoElem>, D, D>: SliceArg<D, OutDim = D>,
{
    unsafe {
        let r = axis_slice_unsafe(a, None, None, Some(-1), axis, a_ndim);
        debug_assert!(r.is_ok());
        r.unwrap_unchecked()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ndarray::{array, Array, ArrayD, IxDyn};

    /// Tests on IxN arrays.
    #[test]
    fn axis_slice_doc() {
        let a = array![[1, 2, 3], [4, 5, 6], [7, 8, 9]];

        assert_eq!(
            axis_slice(&a, Some(0), Some(1), Some(1), Some(1)).unwrap(),
            array![[1], [4], [7]]
        );
        assert_eq!(
            axis_slice(&a, Some(0), Some(2), Some(1), Some(0)).unwrap(),
            array![[1, 2, 3], [4, 5, 6]]
        );
    }

    /// Tests on IxN arrays with negative step.
    #[test]
    fn axis_slice_neg_step() {
        let a = array![[1, 2, 3, 4, 5], [0, 1, 4, 9, 16]];
        assert_eq!(
            axis_slice(&a, Some(2), Some(0), Some(-1), None).unwrap(),
            array![[3, 2], [4, 1]]
        );
        assert_eq!(
            axis_slice(&a, Some(-2), Some(-4), Some(-1), None).unwrap(),
            array![[4, 3], [9, 4]]
        );
    }

    /// Tests on IxN arrays with negative indices.
    #[test]
    fn axis_slice_neg_indices_weird() {
        let a = array![1, 2, 3, 4];
        assert_eq!(
            unsafe { axis_slice(&a, Some(-2), Some(-5), Some(-1), None) }.unwrap(),
            array![3, 2, 1]
        );
        assert_eq!(
            unsafe { axis_slice_unsafe(&a, Some(-2), Some(-5), Some(-1), 0, a.ndim()) }.unwrap(),
            array![3, 2, 1]
        );
    }

    /// Test on IxDyn Arrays.
    #[test]
    fn axis_slice_doc_dyn() {
        let a = {
            let mut y: Array<_, IxDyn> = ArrayD::<i64>::zeros(IxDyn(&[2, 3]));
            y[[0, 0]] = 5;
            y[[0, 1]] = 6;
            y[[0, 2]] = 7;
            y[[1, 0]] = 1;
            y[[1, 1]] = 2;
            y[[1, 2]] = 3;

            y
        };

        assert_eq!(
            axis_slice(&a, Some(0), Some(1), Some(1), Some(1))
                .unwrap()
                .into_dimensionality()
                .unwrap(),
            array![[5], [1]]
        );
    }
}
