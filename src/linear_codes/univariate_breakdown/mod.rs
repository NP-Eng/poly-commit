use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{merkle_tree::Config, sponge::CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;

use super::{BreakdownPCParams, LinearEncode};

mod tests;

/// The univariate Breakdown polynomial commitment scheme based on [[Breakdown]][bd].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [bd]: https://eprint.iacr.org/2021/1043.pdf
pub struct UnivariateBreakdown<
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    H: CRHScheme,
> {
    _phantom: PhantomData<(F, C, S, P, H)>,
}

impl<F, C, S, P, H> LinearEncode<F, C, P, H> for UnivariateBreakdown<F, C, S, P, H>
where
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    P::Point: Into<F>,
    H: CRHScheme,
{
    type LinCodePCParams = BreakdownPCParams<F, C, H>;

    fn setup<R: RngCore>(
        rng: &mut R,
        leaf_hash_params: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
        col_hash_params: H::Parameters,
    ) -> Self::LinCodePCParams {
        Self::LinCodePCParams::new(
            rng,
            128,
            (1, 5),
            (82, 1000),
            (10, 6),
            30,
            1 << 20,
            true,
            leaf_hash_params,
            two_to_one_params,
            col_hash_params,
        )
    }

    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Vec<F> {
        todo!()
    }

    /// For a univariate polynomial, we simply return the list of coefficients
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
}

impl<F, C, S, P, H> UnivariateBreakdown<F, C, S, P, H>
where
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    P::Point: Into<F>,
    H: CRHScheme,
{
}
