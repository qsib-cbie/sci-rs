#[cfg(feature = "unstable")]
use heapless::Vec;
use nalgebra::{Complex, RealField};

use super::Sos;

pub enum FilterOutputType {
    Ba,
    Zpk,
    Sos,
}

#[cfg(feature = "unstable")]
#[derive(Debug)]
pub struct BaFormatFilterSt<F: RealField, const N: usize> {
    pub b: [F; N],
    pub a: [F; N],
}

#[cfg(feature = "unstable")]
#[derive(Debug)]
pub struct ZpkFormatFilterSt<F: RealField + Copy, const N: usize>
where
    Vec<Complex<F>, N>: Sized,
{
    pub z: Vec<Complex<F>, N>,
    pub p: Vec<Complex<F>, N>,
    pub k: F,
}

#[cfg(feature = "unstable")]
#[derive(Debug)]
pub struct SosFormatFilterSt<F: RealField + Copy, const N: usize> {
    pub sos: Vec<Sos<F>, N>,
}

#[cfg(feature = "unstable")]
#[derive(Debug)]
pub enum DigitalFilterSt<F: RealField + Copy + Sized, const N: usize>
where
    [(); { N * 2 + 1 }]: Sized,
    [(); { N * 2 }]: Sized,
{
    Ba(BaFormatFilterSt<F, { N * 2 + 1 }>),
    Zpk(ZpkFormatFilterSt<F, { N * 2 }>),
    Sos(SosFormatFilterSt<F, N>),
}

#[cfg(feature = "unstable")]
impl<F: RealField + Copy, const N: usize> ZpkFormatFilterSt<F, N> {
    pub fn new(z: Vec<Complex<F>, N>, p: Vec<Complex<F>, N>, k: F) -> Self {
        ZpkFormatFilterSt { z, p, k }
    }
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

    #[cfg(feature = "unstable")]
    #[test]
    fn can_create_digital_filter() {
        let single_side = DigitalFilter::<f32, 1>::Sos(SosFormatFilter { sos: Vec::new() });
        match single_side {
            DigitalFilter::Sos(s) => {
                assert_eq!(s.sos.capacity(), 1);
            }
            _ => unreachable!(),
        }

        let single_side = DigitalFilter::<f32, 1>::Zpk(ZpkFormatFilter {
            z: Vec::new(),
            p: Vec::new(),
            k: 1.,
        });
        match single_side {
            DigitalFilter::Zpk(zpk) => {
                assert_eq!(zpk.z.capacity(), 2);
                assert_eq!(zpk.p.capacity(), 2);
            }
            _ => unreachable!(),
        }
    }

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
