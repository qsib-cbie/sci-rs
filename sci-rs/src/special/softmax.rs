use ndarray::{Array, ArrayBase, Axis};
use num_traits::Float;

/// Compute the softmax function along the specified axis.
///
/// The softmax function is defined as:
/// ```text
/// softmax(x_i) = exp(x_i) / sum(exp(x_j) for j in axis)
/// ```
///
/// This function is usually used in machine learning to normalize the output of a neural network to a probability
/// distribution.
/// ```
/// use ndarray::{array, Axis};
/// use sci_rs::special::softmax;
///
/// let a = array![[1., 2., 3.], [4., 5., 6.0_f32]];
/// let b = softmax(&a, Axis(0)).mapv(|x| (x * 100.0).round() / 100.0);
/// assert_eq!(b, array![[0.05, 0.05, 0.05], [0.95, 0.95, 0.95]]);
/// let c = softmax(&a, Axis(1)).mapv(|x| (x * 100.0).round() / 100.0);
/// assert_eq!(c, array![[0.09, 0.24, 0.67], [0.09, 0.24, 0.67]]);
/// ```
///
/// # Arguments
///
/// * `x`: The input array.
/// * `axis`: The axis along which to compute the softmax function (so every slice along the axis will sum to 1).
pub fn softmax<A, S, D>(x: &ArrayBase<S, D>, axis: Axis) -> Array<A, D>
where
    A: Float + 'static,
    S: ndarray::Data<Elem = A>,
    D: ndarray::RemoveAxis,
{
    let mut res = Array::uninit(x.raw_dim());
    for (arr, mut res) in x.lanes(axis).into_iter().zip(res.lanes_mut(axis)) {
        let max = arr
            .iter()
            // If we have NaN and the comparison fails, the max can be arbitrary as the sum and the whole result
            // will be NaN anyway, so we use an arbitrary ordering.
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let max = match max {
            Some(max) => *max,
            None => continue,
        };
        let mut sum = A::zero();
        for (i, x) in res.indexed_iter_mut() {
            let v = (arr[i] - max).exp();
            sum = sum + v;
            x.write(v);
        }
        for x in res.iter_mut() {
            // Safety: we wrote to every single element of the `res` array in the previous loop.
            x.write(*unsafe { x.assume_init_ref() } / sum);
        }
    }
    // Safety: we wrote to every single element of the array.
    unsafe { res.assume_init() }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use ndarray::{array, Axis};

    #[test]
    fn softmax() {
        let a = array![[1., 2., 3.], [4., 5., 6.0_f32]];
        let b = super::softmax(&a, Axis(0)).mapv(|x| (x * 100.0).round() / 100.0);
        assert_eq!(b, array![[0.05, 0.05, 0.05], [0.95, 0.95, 0.95]]);
        let c = super::softmax(&a, Axis(1)).mapv(|x| (x * 100.0).round() / 100.0);
        assert_eq!(c, array![[0.09, 0.24, 0.67], [0.09, 0.24, 0.67]]);

        // examples copied from scipy softmax documentation

        let x = array![[1., 0.5, 0.2, 3.], [1., -1., 7., 3.], [2., 12., 13., 3.]];

        let m = super::softmax(&x, Axis(0));
        let y = array![
            [0.211942, 0.00001013, 0.00000275, 0.333333],
            [0.211942, 0.00000226, 0.00247262, 0.333333],
            [0.576117, 0.999988, 0.997525, 0.333333]
        ];
        assert_relative_eq!(m, y, epsilon = 1e-5);

        let m = super::softmax(&x, Axis(1));
        let y = array![
            [1.05877e-01, 6.42177e-02, 4.75736e-02, 7.82332e-01],
            [2.42746e-03, 3.28521e-04, 9.79307e-01, 1.79366e-02],
            [1.22094e-05, 2.68929e-01, 7.31025e-01, 3.31885e-05]
        ];
        assert_relative_eq!(m, y, epsilon = 1e-5);
    }
}
