/// Type of IIR filter
pub enum FilterType {
    /// Butterworth
    /// <https://en.wikipedia.org/wiki/Butterworth_filter>
    Butterworth,
    /// Chebyshev type I
    /// <https://en.wikipedia.org/wiki/Chebyshev_filter>
    ChebyshevI,
    /// Chebyshev type II
    /// <https://en.wikipedia.org/wiki/Chebyshev_filter>
    ChebyshevII,
    /// Cauer/elliptic
    /// <https://en.wikipedia.org/wiki/Elliptic_filter>
    CauerElliptic,
    /// Bessel/Thomson
    /// <https://en.wikipedia.org/wiki/Bessel_filter>
    BesselThomson(BesselThomsonNorm),
}

/// Bessel-Thomson filter normalization
pub enum BesselThomsonNorm {
    /// Phase
    Phase,
    /// Delay
    Delay,
    /// Magnitude
    Mag,
}

/// Type of IIR or FIR filter
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FilterBandType {
    /// Lowpass
    /// (single-sided)
    Lowpass,
    /// Highpass
    /// (single-sided)
    Highpass,
    /// Bandpass
    /// (double-sided)
    Bandpass,
    /// Bandstop or notch filter
    /// (double-sided)
    Bandstop,
}
