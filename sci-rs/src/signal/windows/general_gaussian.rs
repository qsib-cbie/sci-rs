use super::{extend, len_guard, truncate};
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

/// Collection of arguments for window `GeneralCosine` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralGaussian<F>
where
    F: Real,
{
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Shape parameter.
    pub p: F,
    /// The standard deviation, σ.
    pub sigma: F,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl<F> GeneralGaussian<F>
where
    F: Real,
{
    /// Returns a GeneralGaussian struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `p` : float  
    ///     Shape parameter. p = 1 is identical to `gaussian`, p = 0.5 is
    ///     the same shape as the Laplace distribution.
    /// * `sig` : float  
    ///     The standard deviation, σ.
    /// * `sym`:  
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, p: F, sigma: F, sym: bool) -> Self {
        GeneralGaussian { m, p, sigma, sym }
    }
}

#[cfg(feature = "alloc")]
impl<F, W> GetWindow<W> for GeneralGaussian<F>
where
    F: Real,
    W: Real,
{
    /// Return a window with a generalized gaussian shape.
    ///
    /// # Parameters
    /// `self`: [GeneralGaussian]
    ///
    /// # Returns
    /// `w`: `vec<F>`  
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// # Notes
    /// The generalized Gaussian window is defined as  
    /// $$w(n) = e^{ -\frac{1}{2}\left|\frac{n}{\sigma}\right|^{2p} }$$  
    /// the half-power point is at  
    /// $$(2 \log(2))^{1/(2 p)} \sigma$$
    ///
    /// # Example
    /// ```
    /// use sci_rs::signal::windows::{GeneralGaussian, GetWindow};
    /// let window: Vec<f64> = GeneralGaussian::new(51, 1.5, 7., true).get_window();
    /// ```
    ///
    /// This is equivalent to the Python code:
    /// ```custom,{class=language-python}
    /// from scipy import signal
    /// window = signal.windows.general_gaussian(51, p=1.5, sig=7)
    /// ```
    ///
    /// # References
    /// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.general_cosine.html>
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        if len_guard(self.m) {
            return Vec::<W>::new();
        }
        let (m, needs_trunc) = extend(self.m, self.sym);

        let n = (0..m).map(|v| {
            W::from(v).unwrap() - (W::from(m).unwrap() - W::from(1).unwrap()) / W::from(2).unwrap()
        });
        let sig = W::from(self.sigma).unwrap();
        let two_p = W::from(self.p).unwrap() * W::from(2).unwrap();
        let w = n
            .into_iter()
            .map(|v| {
                (v / sig)
                    .abs()
                    .powf(two_p)
                    .mul(W::from(-0.5).unwrap())
                    .exp()
            })
            .collect();

        truncate(w, needs_trunc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn general_gaussian_case_a() {
        let gc = GeneralGaussian::new(20, 0.8, 4.1, true);
        let expected = vec![
            0.14688425, 0.2008071, 0.2687279, 0.35164027, 0.44932387, 0.55972797, 0.67829536,
            0.79725628, 0.90478285, 0.98289486, 0.98289486, 0.90478285, 0.79725628, 0.67829536,
            0.55972797, 0.44932387, 0.35164027, 0.2687279, 0.2008071, 0.14688425,
        ];

        assert_vec_eq(expected, gc.get_window());
    }

    #[test]
    fn general_gaussian_case_b() {
        let gc = GeneralGaussian::new(20, 0.8, 4.1, false);
        let expected = vec![
            0.12465962, 0.17219048, 0.23293383, 0.30828863, 0.3987132, 0.5031538, 0.61839303,
            0.73835825, 0.85338105, 0.94904249, 1., 0.94904249, 0.85338105, 0.73835825, 0.61839303,
            0.5031538, 0.3987132, 0.30828863, 0.23293383, 0.17219048,
        ];

        assert_vec_eq(expected, gc.get_window());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
