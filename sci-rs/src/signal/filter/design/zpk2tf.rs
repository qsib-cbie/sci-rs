use super::{cplx::sort_cplx_st, BaFormatFilter};
use core::{cmp::Ordering, f64::consts::PI, iter::Sum, ops::Mul};

use heapless::Vec;
use nalgebra::{
    allocator::Allocator, ArrayStorage, Complex, ComplexField, Const, DefaultAllocator, RealField,
    SMatrix, SVector, Storage, Vector, U1,
};
use num_traits::{Float, Zero};

use crate::signal::filter::design::cplx::cplxreal_st;

use super::{FilterBandType, FilterOutputType, FilterType, Sos, SosFormatFilter, ZpkFormatFilter};

///
/// """
/// Return polynomial transfer function representation from zeros and poles
///
/// Parameters
/// ----------
/// z : array_like
///     Zeros of the transfer function.
/// p : array_like
///     Poles of the transfer function.
/// k : float
///     System gain.
///
/// Returns
/// -------
/// b : ndarray
///     Numerator polynomial coefficients.
/// a : ndarray
///     Denominator polynomial coefficients.
///
/// """
///
pub fn zpk2tf_st<C, F, const N: usize>(
    z: &Vec<C, N>,
    p: &Vec<C, N>,
    k: F,
) -> BaFormatFilter<F, { N + 1 }>
where
    C: ComplexField<RealField = F>,
    F: Float + RealField,
    [(); { N + 1 }]: Sized,
{
    // not possible to have length of shape > 1, but handled in zpk2tf python code

    let b = poly_st(z)
        .into_iter()
        .map(|bi| <C as ComplexField>::from_real(k) * bi)
        .collect::<Vec<_, _>>();
    let a = poly_st(p);

    // Use real output if possible.
    let poly2ba = |x: Vec<C, { N + 1 }>, y: &Vec<C, N>| -> Vec<C, { N + 1 }> {
        let mut x = x;
        let mut pos_roots = y
            .iter()
            .cloned()
            .filter(|i| i.clone().imaginary() > <C as ComplexField>::RealField::zero())
            .collect::<Vec<_, N>>();
        let mut neg_roots = y
            .iter()
            .cloned()
            .filter(|i| i.clone().imaginary() < <C as ComplexField>::RealField::zero())
            .map(|i| i.conjugate())
            .collect::<Vec<_, N>>();
        if pos_roots.len() == neg_roots.len() {
            sort_cplx_st(&mut pos_roots);
            sort_cplx_st(&mut neg_roots);
            if pos_roots.into_iter().zip(neg_roots).all(|(p, n)| p == n) {
                x = x
                    .into_iter()
                    .map(|xi| ComplexField::from_real(xi.real()))
                    .collect::<Vec<_, _>>();
            }
        }
        x
    };

    let bv = poly2ba(b, z);
    let av = poly2ba(a, p);

    let mut b = [F::zero(); { N + 1 }];
    bv.into_iter()
        .zip(b.iter_mut())
        .for_each(|(bvi, bi)| *bi = bvi.real());
    let mut a = [F::zero(); { N + 1 }];
    av.into_iter()
        .zip(a.iter_mut())
        .for_each(|(avi, ai)| *ai = avi.real());

    BaFormatFilter { b, a }
}

pub fn poly_st<F, const N: usize>(z: &Vec<F, N>) -> Vec<F, { N + 1 }>
where
    F: ComplexField,
    [(); { N + 1 }]: Sized,
{
    let mut a: Vec<F, { N + 1 }> = Vec::new();
    a.push(F::one());

    const KER: usize = 2;
    for zi in z {
        let mut b = Vec::new();
        b.resize(a.len() + 1, F::zero());
        let k = [F::one(), -zi.clone()];
        for i in 0..a.len() + KER - 1 {
            let u_i = if i > a.len() { i - KER } else { 0 };
            let u_f = i.min(a.len() - 1);
            if u_i == u_f {
                b[i] += a[u_i].clone() * k[i - u_i].clone();
            } else {
                for u in u_i..(u_f + 1) {
                    if i - u < KER {
                        b[i] += a[u].clone() * k[i - u].clone();
                    }
                }
            }
        }
        a = b;
    }

    let mut roots = z.clone();
    sort_cplx_st(&mut roots);
    let mut root_conjs = z
        .iter()
        .map(|zi| zi.clone().conjugate())
        .collect::<Vec<_, N>>();
    sort_cplx_st(&mut root_conjs);
    if roots.into_iter().zip(root_conjs).all(|(a, b)| a == b) {
        a = a
            .into_iter()
            .map(|ai| ComplexField::from_real(ai.real()))
            .collect();
    }

    a
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn matches_scipy_poly() {
        let mut a: Vec<Complex<f64>, 2> = Vec::new();
        for i in [1., 1.] {
            a.push(Complex::new(i, 0.));
        }
        let c = poly_st(&a);
        [1., -2., 1.].iter().zip(c.iter()).for_each(|(e, a)| {
            assert_relative_eq!(*e, a.real());
            assert_relative_eq!(0., a.imaginary());
        });

        a.clear();
        a.push(Complex::new(0.98924866, 0.03710237));
        a.push(Complex::new(0.98924866, -0.03710237));
        let c = poly_st(&a);
        [1., -1.97849731, 0.97998949]
            .iter()
            .zip(c.iter())
            .for_each(|(e, a)| {
                assert_relative_eq!(*e, a.real(), max_relative = 1e-7);
                assert_relative_eq!(0., a.imaginary(), max_relative = 1e-7);
            });

        a.clear();
        a.push(Complex::new(0.96189799, 0.03364097));
        a.push(Complex::new(0.96189799, -0.03364097));
        let c = poly_st(&a);
        [1., -1.92379599, 0.92637947]
            .iter()
            .zip(c.iter())
            .for_each(|(e, a)| {
                assert_relative_eq!(*e, a.real(), max_relative = 1e-7);
                assert_relative_eq!(0., a.imaginary(), max_relative = 1e-7);
            });
    }
}
