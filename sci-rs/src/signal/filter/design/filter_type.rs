pub enum FilterType {
    Butterworth,
    ChebyshevI,
    ChebyshevII,
    CauerElliptic,
    BesselThomson(BesselThomsonNorm),
}

pub enum BesselThomsonNorm {
    Phase,
    Delay,
    Mag,
}

pub enum FilterBandType {
    Lowpass,
    Highpass,
    Bandpass,
    Bandstop,
}
