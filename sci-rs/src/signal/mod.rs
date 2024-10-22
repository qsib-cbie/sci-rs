/// Digital Filtering
pub mod filter;

/// Signal Generation
pub mod wave;

/// Convolution
#[cfg(feature = "std")]
pub mod convolve;

/// Window functions  
/// This contains all window functions in the
/// [scipy.signal.windows](https://docs.scipy.org/doc/scipy/reference/signal.windows.html#module-scipy.signal.windows)
/// namespace.  
/// The convenience function
/// [`get_windows`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window)
/// in the [scipy.signal](https://docs.scipy.org/doc/scipy/reference/signal.html#window-functions)
/// namespace is located here.
pub mod windows;

/// Signal Resampling
#[cfg(feature = "std")]
pub mod resample;
