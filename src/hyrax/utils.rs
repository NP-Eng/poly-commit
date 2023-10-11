
/// Transforms a flat vector into a matrix in column major order.
/// For example, if flat = [1, 2, 3, 4, 5, 6] and n = 2, m = 3, then
/// the output is [[1, 3, 5], [2, 4, 6]].
pub(crate) fn flat_to_matrix_column_major<T>(flat: &[T], n: usize, m: usize) -> Vec<Vec<T>> {
    assert_eq!(flat.len(), n * m, "n * m should coincide with the length of flat");
    let mut res = vec![vec![0; m]; n];
    for row in 0..n {
        for col in 0..m {
            res[row][col] = flat[col * n + row];
        }
    }
    res
}