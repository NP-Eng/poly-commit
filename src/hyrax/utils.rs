use ark_ff::Field;

/// Transforms a flat vector into a matrix in column major order.
/// For example, if flat = [1, 2, 3, 4, 5, 6] and n = 2, m = 3, then
/// the output is [[1, 3, 5], [2, 4, 6]].
pub(crate) fn flat_to_matrix_column_major<T: Copy>(flat: &[T], n: usize, m: usize) -> Vec<Vec<T>> {
    assert_eq!(
        flat.len(),
        n * m,
        "n * m should coincide with the length of flat"
    );
    let mut res = Vec::new();

    for row in 0..n {
        res.push((0..m).map(|col| flat[col * n + row]).collect())
    }
    res
}

// This function computes all evaluations of the MLE EQ(i, values) for i
// between 0...0 and 1...1 (n-bit strings). This results in essentially
// the same as the tensor_inner function in the multilinear_ligero module,
// the difference being the endianness of the order of the output.
pub(crate) fn tensor_prime<F: Field>(values: &[F]) -> Vec<F> {
    if values.len() == 0 {
        return vec![F::one()];
    }

    let tail = tensor_prime(&values[1..]);

    tail.iter()
        .map(|v| *v * (F::one() - values[0]))
        .chain(tail.iter().map(|v| *v * values[0]))
        .collect()
}
