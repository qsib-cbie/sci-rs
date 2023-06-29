use nalgebra::{Complex, RealField};

use super::Sos;

pub enum FilterOutputType {
    Ba,
    Zpk,
    Sos,
}

#[cfg(feature = "use_std")]
#[derive(Debug)]
pub struct BaFormatFilter<F: RealField> {
    pub b: Vec<F>,
    pub a: Vec<F>,
}

#[cfg(feature = "use_std")]
#[derive(Debug)]
pub struct ZpkFormatFilter<F: RealField + Copy> {
    pub z: Vec<Complex<F>>,
    pub p: Vec<Complex<F>>,
    pub k: F,
}

#[cfg(feature = "use_std")]
#[derive(Debug)]
pub struct SosFormatFilter<F: RealField + Copy> {
    pub sos: Vec<Sos<F>>,
}

#[cfg(feature = "use_std")]
#[derive(Debug)]
pub enum DigitalFilter<F: RealField + Copy + Sized> {
    Ba(BaFormatFilter<F>),
    Zpk(ZpkFormatFilter<F>),
    Sos(SosFormatFilter<F>),
}

#[cfg(feature = "use_std")]
impl<F: RealField + Copy> ZpkFormatFilter<F> {
    pub fn new(z: Vec<Complex<F>>, p: Vec<Complex<F>>, k: F) -> Self {
        ZpkFormatFilter { z, p, k }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "use_std")]
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
