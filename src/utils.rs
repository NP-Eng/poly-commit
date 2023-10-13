#[cfg(not(feature = "std"))]
use num_traits::Float;

/// Entropy function
fn _ent(x: f64) -> f64 {
    assert!(0f64 <= x && x <= 1f64);
    if x == 0f64 || x == 1f64 {
        0f64
    } else {
        -x * x.log2() - (1.0 - x) * (1.0 - x).log2()
    }
}

/// ceil of a * b, where a is integer and b is a rational number
#[inline]
fn _ceil_mul(a: usize, b: (usize, usize)) -> usize {
    (a * b.0 + b.1 - 1) / b.1
}

/// Return ceil(x / y).
pub(crate) fn ceil_div(x: usize, y: usize) -> usize {
    // XXX. warning: this expression can overflow.
    (x + y - 1) / y
}
