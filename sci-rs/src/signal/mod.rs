/// Digital Filtering  
/// Contains functions from [Filtering section of
/// `scipy.signal`](https://docs.scipy.org/doc/scipy/reference/signal.html#filtering).
pub mod filter;

/// Signal Generation  
/// Contains functions from the [Waveforms section of
/// `scipy.signal`](<https://docs.scipy.org/doc/scipy/reference/signal.html#waveforms>).
pub mod wave;

/// Convolution  
/// Contains functions from the [Convolution section of
/// `scipy.signal`](<https://docs.scipy.org/doc/scipy/reference/signal.html#convolution>).
#[cfg(feature = "std")]
pub mod convolve;

/// Window functions  
/// This contains all window functions in the
/// [`scipy.signal.windows`](https://docs.scipy.org/doc/scipy/reference/signal.windows.html#module-scipy.signal.windows)
/// namespace.  
/// The convenience function
/// [`get_windows`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window)
/// in the [scipy.signal](https://docs.scipy.org/doc/scipy/reference/signal.html#window-functions)
/// namespace is located here.
pub mod windows;

// Convenience macro for down-stream users
/// This provides the same convenience function as in the scipy.signal namespace. Returns a window
/// of given length and type.
///
/// # Parameters:
/// - `window`: The type of window to create.  
/// - `m`: Number of samples in the window.  
/// - `fftbins`: If true (default), creates a "periodic" window, ready to use with ifftshift and be
///   multiplied by the result of an FFT. If False, create a "symmetric" window, for use in filter design.
#[doc(inline)]
pub use crate::_signal_windows_getWindow as get_window;

/// Signal Resampling  
/// This contains only the
/// [`resample`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample.html#scipy.signal.resample)
/// function from `scipy.signal`.
#[cfg(feature = "std")]
pub mod resample;
