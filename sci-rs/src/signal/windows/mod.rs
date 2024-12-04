use crate::special;
use nalgebra::RealField;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Corresponding window representation for tuple-structs of [Window] variants.
#[cfg(feature = "alloc")]
pub trait GetWindow<W = f64>
where
    W: Real,
{
    /// Returns a window of given length and type.
    ///
    /// # Parameters
    /// `self`:  
    ///     The type of window to construct, typically consists of at least the following
    ///     arguments:
    /// * `Nx`: usize  
    ///     The number of samples in the window
    /// * `fftbins/~sym`: bool  
    ///     If fftbins=true/sym=false, create a “periodic” window, ready to use with ifftshift and
    ///     be multiplied by the result of an FFT. This is the default behaviour in If
    ///     fftbins=false/sym=true, create a "symmetric" window, for use in filter design.
    /// * `*args`:  
    ///     Other arguments relevant to the window type.
    ///     
    /// # Reference
    /// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.get_window.html>
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W>;
}

/// Private function for windows implementing [GetWindow]
/// Handle small or incorrect window lengths.
#[inline(always)]
fn len_guard(m: usize) -> bool {
    m <= 1
}

/// Private function for windows implementing [GetWindow]
/// Extend window by 1 sample if needed for DFT-even symmetry.
#[inline(always)]
fn extend(m: usize, sym: bool) -> (usize, bool) {
    if !sym {
        (m + 1, true)
    } else {
        (m, false)
    }
}

/// Private function for windows implementing [GetWindow]
/// Truncate window by 1 sample if needed for DFT-even symmetry.
#[inline(always)]
fn truncate<W>(mut w: Vec<W>, needed: bool) -> Vec<W> {
    if needed {
        w.pop();
    }
    w
}

mod blackman;
mod boxcar;
mod general_cosine;
mod general_gaussian;
mod general_hamming;
mod hamming;
mod kaiser;
mod nuttall;
mod triangle;
pub use blackman::Blackman;
pub use boxcar::Boxcar;
pub use general_cosine::GeneralCosine;
pub use general_gaussian::GeneralGaussian;
pub use general_hamming::GeneralHamming;
pub use hamming::Hamming;
pub use kaiser::Kaiser;
pub use nuttall::Nuttall;
pub use triangle::Triangle;

/// This collects all structs that implement the [GetWindow] trait.  
/// This allows for running `.get_window()` on the struct, which can then be, for example, used in
/// Firwin.
// Ordering is as in accordance with
// https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html.
#[derive(Debug, Clone, PartialEq)]
// Is it possible for the enums which wraps the structs to only require the generic that the struct
// implements GetWindow?
pub enum Window<F>
where
    F: Real,
{
    /// [Boxcar] window, also known as a rectangular window or Dirichlet window; This is equivalent
    /// to no window at all.
    Boxcar(Boxcar),
    /// [Triangle] window.
    Triangle(Triangle),
    /// [Blackman] window.
    Blackman(Blackman),
    /// [Hamming] window.
    Hamming(Hamming),
    // Hann,
    // Bartlett,
    // Flattop,
    // Parzen,
    // Bohman,
    // BlackmanHarris,
    /// [Nuttall] window.
    Nuttall(Nuttall),
    // BartHann,
    // Cosine,
    // Exponential,
    // Tukey,
    // Taylor,
    // Lanczos,
    /// [Kaiser] window.
    // Needs Beta
    Kaiser(Kaiser<F>),
    // KaiserBesselDerived, // Needs Beta
    // Gaussian, // Needs Standard Deviation
    /// [GeneralCosine] window, a generic weighted sum of cosine term windows.
    // Needs Weighting Coefficients
    GeneralCosine(GeneralCosine<F>),
    /// [GeneralGaussian] window.
    // Needs Power, Width
    GeneralGaussian(GeneralGaussian<F>),
    /// [GeneralHamming] window.
    // Needs Window Coefficients.
    GeneralHamming(GeneralHamming<F>),
    // Dpss, // Needs Normalized Half-Bandwidth.
    // Chebwin, // Needs Attenuation.
}

impl<F, W> GetWindow<W> for Window<F>
where
    F: Real,
    W: Real + Float + RealField + special::Bessel,
{
    fn get_window(&self) -> Vec<W> {
        match &self {
            Window::Boxcar(x) => x.get_window(),
            Window::Triangle(x) => x.get_window(),
            Window::Blackman(x) => x.get_window(),
            Window::Hamming(x) => x.get_window(),
            Window::Nuttall(x) => x.get_window(),
            Window::Kaiser(x) => x.get_window(),
            Window::GeneralCosine(x) => x.get_window(),
            Window::GeneralGaussian(x) => x.get_window(),
            Window::GeneralHamming(x) => x.get_window(),
        }
    }
}

/// This provides a set of enum variants that for use in [get_window].
#[derive(Debug, Clone, PartialEq)] // Derive eq?
pub enum GetWindowBuilder<'a, F>
where
    F: Real,
{
    /// [Boxcar] window, also known as a rectangular window or Dirichlet window; This is equivalent
    /// to no window at all.
    Boxcar,
    /// [Triangle] window.
    Triangle,
    /// [Blackman] window.
    Blackman,
    /// [Hamming] window.
    Hamming,
    // Hann,
    // Bartlett,
    // Flattop,
    // Parzen,
    // Bohman,
    // BlackmanHarris,
    /// [Nuttall] window.
    Nuttall,
    // BartHann,
    // Cosine,
    // Exponential,
    // Tukey,
    // Taylor,
    // Lanczos,
    /// [Kaiser] window.
    Kaiser {
        /// Shape parameter `β`, please refer to [Kaiser].
        beta: F,
    },
    // KaiserBesselDerived, // Needs Beta
    // Gaussian, // Needs Standard Deviation
    /// [GeneralCosine] window: Generic weighted sum of cosine term windows.
    GeneralCosine {
        /// Weighting Coefficients `a`, please refer to [GeneralCosine].
        weights: &'a [F],
    },
    /// [GeneralGaussian] window: Generalized Gaussian Shape.
    GeneralGaussian {
        /// Shape parameter.
        p: F,
        /// The standard deviation, σ.
        width: F,
    },
    /// [GeneralHamming] window.
    // Needs Window Coefficients.
    GeneralHamming {
        /// Window coefficient, ɑ
        coefficient: F,
    },
    // Dpss, // Needs Normalized Half-Bandwidth.
    // Chebwin, // Needs Attenuation.
}

/// Return a window of a given length and type.
///
/// Parameters
/// ----------
/// * `window`: [GetWindowBuilder]  
///     The type of window to create. See below for more details.
/// * `Nx`: usize  
///     The number of samples in the window.
/// * `fftbins`: bool, optional  
///     If True (default), create a "periodic" window, ready to use with `ifftshift` and be
///     multiplied by the result of an FFT (see also :func:`~scipy.fft.fftfreq`).  
///     If False, create a "symmetric" window, for use in filter design.
///
/// Returns
/// -------
/// * `get_window` : ndarray
///     Returns a window of length `Nx` and type `window`
///
/// Notes
/// -----
/// Window types:
/// * [Boxcar]
/// * [Triangle]
/// * [Blackman]
/// * [Hamming]
// Hann,
// Bartlett,
// Flattop,
// Parzen,
// Bohman,
// BlackmanHarris,
/// * [Nuttall]
// BartHann,
// Cosine,
// Exponential,
// Tukey,
// Taylor,
// Lanczos,
/// * [Kaiser] // Needs Beta
// KaiserBesselDerived, // Needs Beta
// Gaussian, // Needs Standard Deviation
/// * [GeneralCosine]
/// * [GeneralGaussian] // Needs Power, Width
/// * [GeneralHamming] // Needs Window Coefficients.
// Dpss, // Needs Normalized Half-Bandwidth.
// Chebwin, // Needs Attenuation.
///
/// Examples
/// -----
/// ```
/// use approx:: assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::{get_window, GetWindow, GetWindowBuilder};
///
/// let window_struct = get_window(GetWindowBuilder::<f64>::Hamming, 3, None);
/// let window: Vec<f64> = window_struct.get_window();
/// let expected = vec![0.08, 0.77, 0.77];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
///
/// # References
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.get_window.html>
pub fn get_window<F>(window: GetWindowBuilder<'_, F>, nx: usize, fftbins: Option<bool>) -> Window<F>
where
    F: Real,
{
    match window {
        GetWindowBuilder::Boxcar => Window::Boxcar(Boxcar {
            m: nx,
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::Triangle => Window::Triangle(Triangle {
            m: nx,
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::Blackman => Window::Blackman(Blackman {
            m: nx,
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::Hamming => Window::Hamming(Hamming {
            m: nx,
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::Nuttall => Window::Nuttall(Nuttall {
            m: nx,
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::Kaiser { beta } => Window::Kaiser(Kaiser {
            m: nx,
            beta,
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::GeneralCosine { weights } => Window::GeneralCosine(GeneralCosine {
            m: nx,
            a: weights.into(),
            sym: !fftbins.unwrap_or(true),
        }),
        GetWindowBuilder::GeneralGaussian { p, width } => {
            Window::GeneralGaussian(GeneralGaussian {
                m: nx,
                p,
                sigma: width,
                sym: !fftbins.unwrap_or(true),
            })
        }
        GetWindowBuilder::GeneralHamming { coefficient } => {
            Window::GeneralHamming(GeneralHamming {
                m: nx,
                alpha: coefficient,
                sym: !fftbins.unwrap_or(true),
            })
        }
    }
}
