use ark_ff::Field;
#[cfg(not(feature = "std"))]
use num_traits::Float;
#[cfg(feature = "parallel")]
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IndexedParallelIterator,
};

/// Entropy function
fn ent(x: f64) -> f64 {
    assert!(0f64 <= x && x <= 1f64);
    if x == 0f64 || x == 1f64 {
        0f64
    } else {
        -x * x.log2() - (1.0 - x) * (1.0 - x).log2()
    }
}

/// ceil of a * b, where a is integer and b is a rational number
#[inline]
fn ceil_mul(a: usize, b: (usize, usize)) -> usize {
    (a * b.0 + b.1 - 1) / b.1
}

/// Return ceil(x / y).
pub(crate) fn ceil_div(x: usize, y: usize) -> usize {
    // XXX. warning: this expression can overflow.
    (x + y - 1) / y
}

#[inline]
pub(crate) fn get_num_bytes(n: usize) -> usize {
    ceil_div((usize::BITS - n.leading_zeros()) as usize, 8)
}

#[inline]
pub(crate) fn inner_product<F: Field>(v1: &[F], v2: &[F]) -> F {
    ark_std::cfg_iter!(v1)
        .zip(v2)
        .map(|(li, ri)| *li * ri)
        .sum()
}
