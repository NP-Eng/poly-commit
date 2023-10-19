use core::borrow::Borrow;

use crate::to_bytes;
use crate::utils::IOPTranscript;
use crate::{utils::ceil_div, Error};

use ark_crypto_primitives::{crh::CRHScheme, merkle_tree::Config};
use ark_ff::{FftField, Field, PrimeField};

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::string::ToString;
use ark_std::vec::Vec;

use digest::Digest;
#[cfg(not(feature = "std"))]
use num_traits::Float;

/// This is CSR format
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct SprsMat<F: Field> {
    pub(crate) n: usize,
    pub(crate) m: usize,
    pub(crate) d: usize,
    ind_ptr: Vec<usize>,
    col_ind: Vec<usize>,
    val: Vec<F>,
}

impl<F: Field> SprsMat<F> {
    /// Calulates v.M
    pub(crate) fn row_mul(&self, v: &[F]) -> Vec<F> {
        (0..self.m)
            .map(|j| {
                let ij = self.ind_ptr[j]..self.ind_ptr[j + 1];
                self.col_ind[ij.clone()]
                    .iter()
                    .zip(&self.val[ij])
                    .map(|(&idx, x)| v[idx] * x)
                    .sum::<F>()
            })
            .collect::<Vec<_>>()
    }
    /// m is nrow, n is ncol, d is NNZ in each column
    pub fn new_from_flat(n: usize, m: usize, d: usize, list: &[F]) -> Self {
        let nnz = d * n;
        let mut ind_ptr = vec![0; m + 1];
        let mut col_ind = Vec::<usize>::with_capacity(nnz);
        let mut val = Vec::<F>::with_capacity(nnz);
        for i in 0..m {
            for (c, &v) in list[i * n..(i + 1) * n].iter().enumerate() {
                if v != F::zero() {
                    ind_ptr[i + 1] += 1;
                    col_ind.push(c);
                    val.push(v);
                }
            }
            ind_ptr[i + 1] += ind_ptr[i]
        }
        assert!(ind_ptr[m] <= nnz);
        Self {
            n,
            m,
            d,
            ind_ptr,
            col_ind,
            val,
        }
    }
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
pub(crate) fn get_num_bytes(n: usize) -> usize {
    ceil_div((usize::BITS - n.leading_zeros()) as usize, 8)
}

#[inline]
pub(crate) fn hash_column<F, C, H>(array: Vec<F>, params: &H::Parameters) -> Result<C::Leaf, Error>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
    Vec<F>: Borrow<<H as CRHScheme>::Input>,
    C::Leaf: Sized,
    H::Output: Into<C::Leaf>,
{
    H::evaluate(params, array)
        .map_err(|_| Error::HashingError)
        .map(|x| x.into())
}

/// Generate `t` (not necessarily distinct) random points in `[0, n)` using the current state of `transcript`
pub(crate) fn get_indices_from_transcript<F: PrimeField>(
    n: usize,
    t: usize,
    transcript: &mut IOPTranscript<F>,
) -> Result<Vec<usize>, Error> {
    let bytes_to_squeeze = get_num_bytes(n);
    let mut indices = Vec::with_capacity(t);
    for _ in 0..t {
        let mut bytes: Vec<u8> = vec![0; bytes_to_squeeze];
        transcript
            .get_and_append_byte_challenge(b"i", &mut bytes)
            .map_err(|_| Error::TranscriptError)?;

        // get the usize from Vec<u8>:
        let ind = bytes.iter().fold(0, |acc, &x| (acc << 8) + x as usize);
        // modulo the number of columns in the encoded matrix
        indices.push(ind % n);
    }
    Ok(indices)
}

#[inline]
pub(crate) fn calculate_t<F: PrimeField>(
    sec_param: usize,
    distance: (usize, usize),
    codeword_len: usize,
) -> Result<usize, Error> {
    // Took from the analysis by BCI+20 and Ligero
    // We will find the smallest $t$ such that
    // $(1-\delta)^t + (\rho+\delta)^t + \frac{n}{F} < 2^{-\lambda}$.
    // With $\delta = \frac{1-\rho}{2}$, the expreesion is
    // $2 * (\frac{1+\rho}{2})^t + \frac{n}{F} < 2^(-\lambda)$.

    let field_bits = F::MODULUS_BIT_SIZE as i32;
    let sec_param = sec_param as i32;

    let residual = codeword_len as f64 / 2.0_f64.powi(field_bits);
    let rhs = (2.0_f64.powi(-sec_param) - residual).log2();
    if !(rhs.is_normal()) {
        return Err(Error::InvalidParameters("For the given codeword length and the required security guarantee, the field is not big enough.".to_string()));
    }
    let nom = rhs - 1.0;
    let denom = (1.0 - 0.5 * distance.0 as f64 / distance.1 as f64).log2();
    if !(denom.is_normal()) {
        return Err(Error::InvalidParameters(
            "The distance is wrong".to_string(),
        ));
    }
    let t = (nom / denom).ceil() as usize;
    Ok(if t < codeword_len { t } else { codeword_len }) // This is the `t`
}

/// Only needed for benches and tests
pub struct LeafIdentityHasher;

impl CRHScheme for LeafIdentityHasher {
    type Input = Vec<u8>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        Ok(input.borrow().to_vec().into())
    }
}

/// Only needed for benches and tests
pub struct FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    _phantom: PhantomData<(F, D)>,
}

impl<F, D> CRHScheme for FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    type Input = Vec<F>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_rng: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _parameters: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        let mut dig = D::new();
        dig.update(to_bytes!(input.borrow()).unwrap());
        Ok(dig.finalize().to_vec())
    }
}

#[cfg(test)]
pub(crate) mod tests {

    use crate::utils::to_field;

    use super::*;
    use ark_bls12_377::Fr;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        Polynomial,
    };
    use ark_std::test_rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    // Define some shared testing hashers for univariate & multilinear ligero.

    #[test]
    fn test_sprs_row_mul() {
        let mat: Vec<Fr> = to_field(vec![10, 23, 55, 100, 1, 58, 4, 0, 9]);

        let mat = SprsMat::new_from_flat(3, 3, 3, &mat);
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
