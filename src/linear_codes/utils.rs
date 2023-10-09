use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
use ark_ff::{FftField, Field, PrimeField};

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::vec::Vec;
use ark_std::{string::ToString, test_rng};
#[cfg(not(feature = "std"))]
use num_traits::Float;
#[cfg(feature = "parallel")]
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IndexedParallelIterator,
};

use crate::streaming_kzg::ceil_div;
use crate::Error;

#[derive(Debug)]
pub struct Matrix<F: Field> {
    pub(crate) n: usize,
    pub(crate) m: usize,
    entries: Vec<Vec<F>>,
}

impl<F: Field> Matrix<F> {
    /// Returns a Matrix of dimensions n x m given a list of n * m field elements.
    /// The list should be ordered row-first, i.e. [a11, ..., a1m, a21, ..., a2m, ...].
    ///
    /// # Panics
    /// Panics if the dimensions do not match the length of the list
    pub(crate) fn new_from_flat(n: usize, m: usize, entry_list: &[F]) -> Self {
        assert_eq!(
            entry_list.len(),
            n * m,
            "Invalid matrix construction: dimensions are {} x {} but entry vector has {} entries",
            n,
            m,
            entry_list.len()
        );

        // TODO more efficient to run linearly?
        let entries: Vec<Vec<F>> = (0..n)
            .map(|row| (0..m).map(|col| entry_list[m * row + col]).collect())
            .collect();

        Self { n, m, entries }
    }

    /// Returns a Matrix given a list of its rows, each in turn represented as a list of field elements.
    ///
    /// # Panics
    /// Panics if the sub-lists do not all have the same length.
    pub(crate) fn new_from_rows(row_list: Vec<Vec<F>>) -> Self {
        let m = row_list[0].len();

        for row in row_list.iter().skip(1) {
            assert_eq!(
                row.len(),
                m,
                "Invalid matrix construction: not all rows have the same length"
            );
        }

        Self {
            n: row_list.len(),
            m,
            entries: row_list,
        }
    }

    /// Returns the entry in position (i, j). **Indexing starts at 0 in both coordinates**,
    /// i.e. the first element is in position (0, 0) and the last one in (n - 1, j - 1),
    /// where n and m are the number of rows and columns, respectively.
    ///
    /// Index bound checks are waived for efficiency and behaviour under invalid indexing is undefined
    #[cfg(test)]
    pub(crate) fn entry(&self, i: usize, j: usize) -> F {
        self.entries[i][j]
    }

    /// Returns self as a list of rows
    pub(crate) fn rows(&self) -> Vec<Vec<F>> {
        self.entries.clone()
    }

    /// Returns self as a list of columns
    pub(crate) fn cols(&self) -> Vec<Vec<F>> {
        (0..self.m)
            .map(|col| (0..self.n).map(|row| self.entries[row][col]).collect())
            .collect()
    }

    /// Returns the product v * self, where v is interpreted as a row vector. In other words,
    /// it returns a linear combination of the rows of self with coefficients given by v.
    ///
    /// Panics if the length of v is different from the number of rows of self.
    pub(crate) fn row_mul(&self, v: &[F]) -> Vec<F> {
        assert_eq!(
            v.len(),
            self.n,
            "Invalid row multiplication: vectir has {} elements whereas each matrix column has {}",
            v.len(),
            self.n
        );

        (0..self.m)
            .map(|col| {
                inner_product(
                    v,
                    &(0..self.n)
                        .map(|row| self.entries[row][col])
                        .collect::<Vec<F>>(),
                )
            })
            .collect()
    }
}

/// Compute the dimensions of an FFT-friendly (over F) matrix with at least n entries.
/// The return pair (n, m) corresponds to the dimensions n x m.
pub(crate) fn compute_dimensions<F: FftField>(n: usize) -> (usize, usize) {
    assert_eq!(
        (n as f64) as usize,
        n,
        "n cannot be converted to f64: aborting"
    );

    let aux = (n as f64).sqrt().ceil() as usize;
    let n_cols = GeneralEvaluationDomain::<F>::new(aux)
        .expect("Field F does not admit FFT with m elements")
        .size();

    (ceil_div(n, n_cols), n_cols)
}

/// Apply reed-solomon encoding to msg.
/// Assumes msg.len() is equal to the order of an FFT domain in F.
/// Returns a vector of length equal to the smallest FFT domain of size at least msg.len() * RHO_INV.
pub(crate) fn reed_solomon<F: FftField>(
    // msg, of length m, is interpreted as a vector of coefficients of a polynomial of degree m - 1
    msg: &[F],
    rho_inv: usize,
) -> Vec<F> {
    let m = msg.len();

    let extended_domain = GeneralEvaluationDomain::<F>::new(m * rho_inv).unwrap_or_else(|| {
        panic!(
            "The field F cannot accomodate FFT for msg.len() * RHO_INV = {} elements (too many)",
            m * rho_inv
        )
    });

    extended_domain.fft(msg)
}

#[inline]
pub(crate) fn inner_product<F: Field>(v1: &[F], v2: &[F]) -> F {
    ark_std::cfg_iter!(v1)
        .zip(v2)
        .map(|(li, ri)| *li * ri)
        .sum()
}

#[inline]
#[cfg(test)]
pub(crate) fn to_field<F: Field>(v: Vec<u64>) -> Vec<F> {
    v.iter().map(|x| F::from(*x)).collect::<Vec<F>>()
}

#[inline]
pub(crate) fn get_num_bytes(n: usize) -> usize {
    ceil_div((usize::BITS - n.leading_zeros()) as usize, 8)
}

// TODO: replace by https://github.com/arkworks-rs/crypto-primitives/issues/112.
#[cfg(test)]
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;

pub(crate) fn poseidon_parameters_for_test<F: PrimeField>() -> PoseidonConfig<F> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut ark = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();

        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        ark.push(res);
    }
    PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, ark, 2, 1)
}

use super::transcript::IOPTranscript;
#[cfg(test)]
pub(crate) fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    use ark_crypto_primitives::sponge::CryptographicSponge;

    PoseidonSponge::new(&poseidon_parameters_for_test())
}

/// Takes as input a struct, and converts them to a series of bytes. All traits
/// that implement `CanonicalSerialize` can be automatically converted to bytes
/// in this manner.
/// From jellyfish lib
#[macro_export]
macro_rules! to_bytes {
    ($x:expr) => {{
        let mut buf = ark_std::vec![];
        ark_serialize::CanonicalSerialize::serialize_compressed($x, &mut buf).map(|_| buf)
    }};
}

/// Generate `t` (not necessarily distinct) random points in `[0, n)` using the current state of `transcript`
pub(crate) fn get_indices_from_transcript<F: PrimeField>(
    n: usize,
    t: usize,
    transcript: &mut IOPTranscript<F>,
) -> Result<Vec<usize>, Error> {
    let mut indices = Vec::with_capacity(t);
    for _ in 0..t {
        let bits = transcript
            .get_and_append_byte_challenge(b"i", n)
            .map_err(|_| Error::TranscriptError)?;

        // get the usize from Vec<u8>:
        let ind = bits.iter().fold(0, |acc, &x| (acc << 1) + x as usize);
        println!("ind: {}", ind);
        indices.push(ind);
    }
    Ok(indices)
}

#[inline]
pub(crate) fn calculate_t<F: PrimeField>(
    sec_param: usize,
    rho_inv: usize,
    codeword_len: usize,
) -> Result<usize, Error> {
    // Took from the analysis by BCI+20 and Ligero
    // We will find the smallest $t$ such that
    // $(1-\delta)^t + (\rho+\delta)^t + \frac{n}{F} < 2^{-\lambda}$.
    // With $\delta = \frac{1-\rho}{2}$, the expreesion is
    // $2 * (\frac{1+\rho}{2})^t + \frac{n}{F} < 2^(-\lambda)$.

    let codeword_len = codeword_len as f64;
    let field_bits = F::MODULUS_BIT_SIZE as i32;
    let sec_param = sec_param as i32;

    let residual = codeword_len / 2.0_f64.powi(field_bits);
    let rhs = (2.0_f64.powi(-sec_param) - residual).log2();
    if !(rhs.is_normal()) {
        return Err(Error::InvalidParameters("For the given codeword length and the required security guarantee, the field is not big enough.".to_string()));
    }
    let nom = rhs - 1.0;
    let denom = (0.5 + 0.5 / rho_inv as f64).log2();
    Ok((nom / denom).ceil() as usize) // This is the `t`
}

#[cfg(test)]
mod tests {

    use ark_bls12_377::Fr;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        Polynomial,
    };
    use ark_std::test_rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    use super::*;

    #[test]
    fn test_matrix_constructor_flat() {
        let entries: Vec<Fr> = to_field(vec![10, 100, 4, 67, 44, 50]);
        let mat = Matrix::new_from_flat(2, 3, &entries);
        assert_eq!(mat.entry(1, 2), Fr::from(50));
    }

    #[test]
    fn test_matrix_constructor_flat_square() {
        let entries: Vec<Fr> = to_field(vec![10, 100, 4, 67]);
        let mat = Matrix::new_from_flat(2, 2, &entries);
        assert_eq!(mat.entry(1, 1), Fr::from(67));
    }

    #[test]
    #[should_panic(expected = "dimensions are 2 x 3 but entry vector has 5 entries")]
    fn test_matrix_constructor_flat_panic() {
        let entries: Vec<Fr> = to_field(vec![10, 100, 4, 67, 44]);
        Matrix::new_from_flat(2, 3, &entries);
    }

    #[test]
    fn test_matrix_constructor_rows() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];
        let mat = Matrix::new_from_rows(rows);
        assert_eq!(mat.entry(2, 0), Fr::from(55));
    }

    #[test]
    #[should_panic(expected = "not all rows have the same length")]
    fn test_matrix_constructor_rows_panic() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58]),
        ];
        Matrix::new_from_rows(rows);
    }

    #[test]
    fn test_cols() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![4, 76]),
            to_field(vec![14, 92]),
            to_field(vec![17, 89]),
        ];

        let mat = Matrix::new_from_rows(rows);

        assert_eq!(mat.cols()[1], to_field(vec![76, 92, 89]));
    }

    #[test]
    fn test_row_mul() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];

        let mat = Matrix::new_from_rows(rows);
        let v: Vec<Fr> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fr
        // we ensure the test will still pass even if Fr changes
        assert_eq!(mat.row_mul(&v), to_field::<Fr>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_encoding() {
        // we use this polynomial to generate the the values we will ask the fft to interpolate

        let rho_inv = 3;
        // `i` is the min number of evaluations we need to interpolate a poly of degree `i - 1`
        for i in 1..10 {
            let deg = (1 << i) - 1;

            let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
            let mut pol = DensePolynomial::rand(deg, rand_chacha);

            while pol.degree() != deg {
                pol = DensePolynomial::rand(deg, rand_chacha);
            }

            let coeffs = &pol.coeffs;

            // size of evals might be larger than deg + 1 (the min. number of evals needed to interpolate): we could still do R-S encoding on smaller evals, but the resulting polynomial will differ, so for this test to work we should pass it in full
            let m = deg + 1;

            let encoded = reed_solomon(&coeffs, rho_inv);

            let large_domain = GeneralEvaluationDomain::<Fr>::new(m * rho_inv).unwrap();

            // the encoded elements should agree with the evaluations of the polynomial in the larger domain
            for j in 0..(rho_inv * m) {
                assert_eq!(pol.evaluate(&large_domain.element(j)), encoded[j]);
            }
        }
    }

    #[test]
    fn test_get_num_bytes() {
        assert_eq!(get_num_bytes(0), 0);
        assert_eq!(get_num_bytes(1), 1);
        assert_eq!(get_num_bytes(9), 1);
        assert_eq!(get_num_bytes(1 << 11), 2);
        assert_eq!(get_num_bytes(1 << 32 - 1), 4);
        assert_eq!(get_num_bytes(1 << 32), 5);
        assert_eq!(get_num_bytes(1 << 32 + 1), 5);
    }
}
