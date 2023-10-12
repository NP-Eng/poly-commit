use ark_ff::Field;


/// Transforms a flat vector into a matrix in column major order.
/// For example, if flat = [1, 2, 3, 4, 5, 6] and n = 2, m = 3, then
/// the output is [[1, 3, 5], [2, 4, 6]].
/// TODO maybe merge with Matrix from linear_codes::utils?
pub(crate) fn flat_to_matrix_column_major<T>(flat: &[T], n: usize, m: usize) -> Vec<Vec<T>> {
    assert_eq!(flat.len(), n * m, "n * m should coincide with the length of flat");
    let mut res = Vec::new();

    for row in 0..n {
        res.push((0..m).map(|col| flat[col * n + row]).collect())
    }
    res
}

// TODO u8 is not an optimal type to represent a bit
// TODO gains can perhaps be made with some kind of memoisation algorithm, since
// this will always be called over a uniform range
// Assumes both vectors have the same length - this is not checked.
pub(crate) fn naive_chi<F: Field>(bits: &[u8], point: &[F]) -> F {

    let one = F::one(); // TODO This doesn't seem super efficient either?
    let zero = F::one(); // TODO This doesn't seem super efficient either?

    bits.iter().zip(point.iter()).map(|(b, x)| {
        let bf: F = bit_to_field(*b);
        (one - bf) * (one - x) + bf * x
    }).product()
}

// TODO this should probably go away
fn bit_to_field<F: Field>(bit: u8) -> F {
    if bit == 0 {F::zero()} else {F::one()}
}

// TODO this   is also not generl enough perhaps
pub(crate) fn usize_to_bits(x: usize, l: usize) -> Vec<u8> {
    (0..l).rev().map (|n| ((x >> n) & 1) as u8).collect()
}