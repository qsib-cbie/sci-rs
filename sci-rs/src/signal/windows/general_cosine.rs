use super::{extend, len_guard, truncate};
use nalgebra::RealField;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

/// Collection of arguments for window `GeneralCosine` for use in [GetWindow].
#[cfg(feature = "alloc")]
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralCosine<F>
where
    F: Real,
{
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Sequence of weighting coefficients.
    pub a: Vec<F>, // Is there a better type, such as impl Into?
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl<F> GeneralCosine<F>
where
    F: Real,
{
    /// Returns a GeneralCosine struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `a`: `Vec<F>`  
    ///     Sequence of weighting coefficients. This uses the convention of being centered on the
    ///     origin, so these will typically all be positive numbers, not alternating sign.
    /// * `sym`:  
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, a: Vec<F>, sym: bool) -> Self {
        GeneralCosine { m, a, sym }
    }
}

#[cfg(feature = "alloc")]
impl<F, W> GetWindow<W> for GeneralCosine<F>
where
    F: Real,
    W: Real + Float + RealField,
{
    /// Return a window of type: GeneralCosine.
    ///
    /// The general cosine window is a generic weighted sum of cosine terms.
    ///
    /// # Parameters
    /// `self`: [GeneralCosine]
    ///
    /// # Returns
    /// `w`: `vec<F>`  
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// # Example
    /// We can create a `flat-top` window named "HFT90D" as following:
    /// ```
    /// use sci_rs::signal::windows::{GeneralCosine, GetWindow};
    ///
    /// let hfd90 = [1., 1.942604, 1.340318, 0.440811, 0.043097].into();
    /// let window: Vec<f64> = GeneralCosine::new(30, hfd90, true).get_window();
    /// ```
    ///
    /// This is equivalent to the Python code:
    /// ```custom,{class=language-python}
    /// from scipy.signal.windows import general_cosine
    /// HFT90D = [1, 1.942604, 1.340318, 0.440811, 0.043097]
    /// window = general_cosine(30, HFT90D, sym=False)
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

        let linspace = (0..m).map(|i| W::from(i).unwrap());
        let fac = linspace.map(|i| W::two_pi() * i / W::from(m - 1).unwrap() - W::pi());
        let w: Vec<_> = self
            .a
            .iter()
            .enumerate()
            .map(|(k, a)| {
                fac.clone()
                    .map(move |f| Float::cos(f * W::from(k).unwrap()) * W::from(*a).unwrap())
            })
            .fold(
                vec![W::from(0).unwrap(); fac.clone().count()],
                // Would this have made more sense if using ndarray::Array1? No, Array1 has no pop.
                |acc, x| acc.iter().zip(x).map(|(&a, b)| a + b).collect(),
            )
            .into_iter()
            .collect();

        truncate(w, needs_trunc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn general_cosine_scipy_eg() {
        // Created with
        // >>> from scipy.signal.windows import general_cosine
        // >>> HFT90D = [1, 1.942604, 1.340318, 0.440811, 0.043097]
        // >>> window = general_cosine(30, HFT90D, sym=False)
        let expected = vec![
            -1.04083409e-16,
            -3.49808964e-03,
            -1.85322176e-02,
            -5.60667246e-02,
            -1.25488810e-01,
            -2.22198500e-01,
            -3.14696393e-01,
            -3.38497088e-01,
            -2.04818447e-01,
            1.72651725e-01,
            8.38783500e-01,
            1.76097559e+00,
            2.81469639e+00,
            3.80321808e+00,
            4.51005597e+00,
            4.76683000e+00,
            4.51005597e+00,
            3.80321808e+00,
            2.81469639e+00,
            1.76097559e+00,
            8.38783500e-01,
            1.72651725e-01,
            -2.04818447e-01,
            -3.38497088e-01,
            -3.14696393e-01,
            -2.22198500e-01,
            -1.25488810e-01,
            -5.60667246e-02,
            -1.85322176e-02,
            -3.49808964e-03,
        ];

        let hfd90 = [1., 1.942604, 1.340318, 0.440811, 0.043097].into();
        let gc = GeneralCosine::new(30, hfd90, false);
        assert_vec_eq(expected, gc.get_window());
    }

    #[test]
    fn constant() {
        let x = 0.42;
        let n = 6;
        let expected = vec![x; n];

        let actual: Vec<f64> = GeneralCosine::new(n, vec![x], true).get_window();
        assert_eq!(expected, actual);
    }

    #[test]
    fn case_b() {
        // Created with
        // >>> from scipy.signal.windows import general_cosine
        // >>> n = 20
        // >>> a = [0.42, 0.50]
        // >>> window = general_cosine(30, a, sym=False)
        let n = 20;
        let a = vec![0.42, 0.50];
        let expected = vec![
            -0.08,
            -0.05552826,
            0.0154915,
            0.12610737,
            0.2654915,
            0.42,
            0.5745085,
            0.71389263,
            0.8245085,
            0.89552826,
            0.92,
            0.89552826,
            0.8245085,
            0.71389263,
            0.5745085,
            0.42,
            0.2654915,
            0.12610737,
            0.0154915,
            -0.05552826,
        ];

        assert_vec_eq(expected, GeneralCosine::new(n, a, false).get_window());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
