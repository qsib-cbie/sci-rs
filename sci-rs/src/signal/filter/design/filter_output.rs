use heapless::Vec;
use nalgebra::{Complex, RealField};

use super::Sos;

pub enum FilterOutputType {
    Ba,
    Zpk,
    Sos,
}

#[derive(Debug)]
pub struct BaFormatFilter<F: RealField, const N: usize> {
    pub b: [F; N],
    pub a: [F; N],
}

#[derive(Debug)]
pub struct ZpkFormatFilter<F: RealField + Copy, const N: usize>
where
    Vec<Complex<F>, N>: Sized,
{
    pub z: Vec<Complex<F>, N>,
    pub p: Vec<Complex<F>, N>,
    pub k: F,
}

#[derive(Debug)]
pub struct SosFormatFilter<F: RealField + Copy, const N: usize> {
    pub sos: Vec<Sos<F>, N>,
}

#[derive(Debug)]
pub enum DigitalFilter<F: RealField + Copy + Sized, const N: usize>
where
    [(); { N * 2 + 1 }]: Sized,
    [(); { N * 2 }]: Sized,
{
    Ba(BaFormatFilter<F, { N * 2 + 1 }>),
    Zpk(ZpkFormatFilter<F, { N * 2 }>),
    Sos(SosFormatFilter<F, N>),
}

impl<F: RealField + Copy, const N: usize> ZpkFormatFilter<F, N> {
    pub fn new(z: Vec<Complex<F>, N>, p: Vec<Complex<F>, N>, k: F) -> Self {
        ZpkFormatFilter { z, p, k }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
}

// pub trait FilterFormat {
//     const ORDER: usize;

//     type Output;

//     fn into(self) -> Self::Output;
// }

// impl<F: RealField + Copy, const N: usize> FilterFormat for BaFormatFilter<F, N> {
//     type Output = Self;

//     const ORDER: usize = N;

//     fn into(self) -> Self::Output {
//         todo!()
//     }
// }

// // pub struct ZpkFormatFilter<F: RealField + Copy + Sized, const N: usize> {
// //     pub z: Vec<Complex<F>, N>,
// //     pub p: Vec<Complex<F>, N>,
// //     pub k: F,
// // }

// // pub struct SosFormatFilter<F: RealField + Copy + Sized, const N: usize> {
// //     pub sos: Vec<Sos<F>, N>,
// // }
