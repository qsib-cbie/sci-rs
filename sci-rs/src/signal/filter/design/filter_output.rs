use nalgebra::{Complex, RealField};

use super::Sos;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Digital filter representation choices
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FilterOutputType {
    /// b/a
    Ba,
    /// zeros/poles/gain
    Zpk,
    /// second order sections
    Sos,
}

/// Numerator/Denominator (b/a) representation of a digital filter
#[cfg(feature = "alloc")]
#[derive(Debug)]
pub struct BaFormatFilter<F: RealField> {
    /// Numerator coefficients
    pub b: Vec<F>,
    /// Denominator coefficients
    pub a: Vec<F>,
}

/// Zeros/Poles/Gain (z/p/k) representation of a digital filter
#[cfg(feature = "alloc")]
#[derive(Debug)]
pub struct ZpkFormatFilter<F: RealField + Copy> {
    /// Zeros of the transfer function.
    pub z: Vec<Complex<F>>,
    /// Poles of the transfer function.
    pub p: Vec<Complex<F>>,
    /// System gain.
    pub k: F,
}

/// Second Order Section Representation
///
/// CMSIS DSP uses negated a1, a2 coefficients as compared to SciPy
/// <https://www.keil.com/pack/doc/CMSIS/DSP/html/group__BiquadCascadeDF1.html>
#[cfg(feature = "alloc")]
#[derive(Debug)]
pub struct SosFormatFilter<F: RealField + Copy> {
    /// Cascaded Second Order Sections
    pub sos: Vec<Sos<F>>,
}

/// Digital Filter Representationin various formats
#[cfg(feature = "alloc")]
#[derive(Debug)]
pub enum DigitalFilter<F: RealField + Copy + Sized> {
    /// Numerator/Denominator (b/a) representation of a digital filter
    Ba(BaFormatFilter<F>),
    /// Zeros/Poles/Gain (z/p/k) representation of a digital filter
    Zpk(ZpkFormatFilter<F>),
    /// Biquad, Second Order Sections (sos) representation of a digital filter
    Sos(SosFormatFilter<F>),
}

#[cfg(feature = "alloc")]
impl<F: RealField + Copy> ZpkFormatFilter<F> {
    /// Create a ZpkFormatFilter
    ///
    /// # Arguments
    /// * `z` - Zeros of the transfer function.
    /// * `p` - Poles of the transfer function.
    /// * `k` - System gain.
    pub fn new(z: Vec<Complex<F>>, p: Vec<Complex<F>>, k: F) -> Self {
        ZpkFormatFilter { z, p, k }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "alloc")]
    #[test]
    fn can_create_digital_filter() {
        let single_side = DigitalFilter::<f32>::Sos(SosFormatFilter { sos: Vec::new() });
        match single_side {
            DigitalFilter::Sos(s) => (),
            _ => unreachable!(),
        }

        let single_side = DigitalFilter::<f32>::Zpk(ZpkFormatFilter {
            z: Vec::new(),
            p: Vec::new(),
            k: 1.,
        });
        match single_side {
            DigitalFilter::Zpk(zpk) => (),
            _ => unreachable!(),
        }
    }
}
