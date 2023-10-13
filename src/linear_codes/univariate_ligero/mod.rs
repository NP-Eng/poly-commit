use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{merkle_tree::Config, sponge::CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_std::borrow::Borrow;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;

use digest::Digest;

use crate::utils::ceil_div;

use super::utils::reed_solomon;
use super::{LigeroPCParams, LinCodeInfo, LinearEncode};

mod tests;

/// The univariate Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
pub struct UnivariateLigero<
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
> {
    _phantom: PhantomData<(F, C, D, S, P)>,
}

impl<F, C, D, S, P> LinearEncode<F, C, D, P> for UnivariateLigero<F, C, D, S, P>
where
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    Vec<u8>: Borrow<C::Leaf>,
    P::Point: Into<F>,
{
    type LinCodePCParams = LigeroPCParams<F, C>;

    fn setup<R>(
        _rng: &mut R,
        leaf_hash_params: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> Self::LinCodePCParams {
        Self::LinCodePCParams::new(128, 4, true, leaf_hash_params, two_to_one_params)
    }

    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Vec<F> {
        reed_solomon(msg, param.rho_inv().0)
    }

    /// For a univariate polynomial, we simply return the list of coefficients.ś
    fn poly_repr(polynomial: &P) -> Vec<F> {
        polynomial.coeffs().to_vec()
    }

    fn point_to_vec(point: P::Point) -> Vec<F> {
        vec![point]
    }

    /// Compute out = [1, z, z^2, ..., z^(n_cols_1)]
    fn tensor(z: &F, left: usize, right: usize) -> (Vec<F>, Vec<F>) {
        let mut left_out = Vec::with_capacity(left);
        let mut pow_a = F::one();
        for _ in 0..left {
            left_out.push(pow_a);
            pow_a *= z;
        }

        let mut right_out = Vec::with_capacity(right);
        let mut pow_b = F::one();
        for _ in 0..right {
            right_out.push(pow_b);
            pow_b *= pow_a;
        }

        (left_out, right_out)
    }

    fn compute_dimensions(n: usize) -> (usize, usize) {
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
}
