use heapless::Vec;
use nalgebra::{Complex, RealField};

use super::Sos;

pub enum FilterOutputType {
    Ba,
    Zpk,
    Sos,
}

pub struct Zpk<F: RealField, const N: usize> {
    pub z: Vec<Complex<F>, N>,
    pub p: Vec<Complex<F>, N>,
    pub k: F,
}

impl<F: RealField, const N: usize> Zpk<F, N> {
    pub fn new(z: Vec<Complex<F>, N>, p: Vec<Complex<F>, N>, k: F) -> Self {
        Zpk { z, p, k }
    }
}

pub enum FilterOutput<F: RealField + Copy + Sized, const N: usize>
where
    [Sos<F>; N / 2 - 1]: Sized,
{
    Ba(([f64; N], [f64; N])),
    Zpk(Zpk<F, N>),
    Sos([Sos<F>; N / 2 - 1]),
}
