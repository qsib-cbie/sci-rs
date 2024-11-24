#[cfg(feature = "alloc")]
use alloc::vec::Vec;
use nalgebra::RealField;
use num_traits::{real::Real, Float};

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
mod triangle;
pub use blackman::Blackman;
pub use boxcar::Boxcar;
pub use general_cosine::GeneralCosine;
pub use general_gaussian::GeneralGaussian;
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
    // Hamming,
    // Hann,
    // Bartlett,
    // Flattop,
    // Parzen,
    // Bohman,
    // BlackmanHarris,
    // Nuttall,
    // BartHann,
    // Cosine,
    // Exponential,
    // Tukey,
    // Taylor,
    // Lanczos,
    // Kaiser Window.
    // Kaiser, // Needs Beta
    // KaiserBesselDerived, // Needs Beta
    // Gaussian, // Needs Standard Deviation
    /// [GeneralCosine] window, a generic weighted sum of cosine term windows.
    // Needs Weighting Coefficients
    GeneralCosine(GeneralCosine<F>),
    /// [GeneralGaussian] window.
    // Needs Power, Width
    GeneralGaussian(GeneralGaussian<F>),
    // GeneralHamming, // Needs Window Coefficients.
    // Dpss, // Needs Normalized Half-Bandwidth.
    // Chebwin, // Needs Attenuation.
}
