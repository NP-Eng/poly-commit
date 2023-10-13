use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    merkle_tree::Config,
    sponge::CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, MultilinearExtension, Polynomial};
use ark_std::borrow::Borrow;
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;
#[cfg(not(feature = "std"))]
use num_traits::Float;

use digest::Digest;

use crate::utils::ceil_div;

use super::{utils::reed_solomon, LigeroPCParams, LinCodeInfo, LinearEncode};

mod tests;

/// The multilinear Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
pub struct MultilinearLigero<
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
> {
    _phantom: PhantomData<(F, C, D, S, P)>,
}

impl<F, C, D, S, P> LinearEncode<F, C, D, P> for MultilinearLigero<F, C, D, S, P>
where
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    Vec<u8>: Borrow<C::Leaf>,
    <P as Polynomial<F>>::Point: Into<Vec<F>>,
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

    fn poly_repr(polynomial: &P) -> Vec<F> {
        polynomial.to_evaluations()
    }

    fn point_to_vec(point: <P as Polynomial<F>>::Point) -> Vec<F> {
        point.into()
    }

    fn tensor(
        point: &<P as Polynomial<F>>::Point,
        left_len: usize,
        _right_len: usize,
    ) -> (Vec<F>, Vec<F>) {
        let point: Vec<F> = Self::point_to_vec(point.clone());

        let split = log2(left_len) as usize;
        let left = &point[..split];
        let right = &point[split..];
        (tensor_inner(left), tensor_inner(right))
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

fn tensor_inner<F: PrimeField>(values: &[F]) -> Vec<F> {
    let one = F::one();
    let anti_values: Vec<F> = values.iter().map(|v| one - *v).collect();

    let mut layer: Vec<F> = vec![one];

    for i in 0..values.len() {
        let mut new_layer = Vec::new();
        for v in &layer {
            new_layer.push(*v * anti_values[i]);
        }
        for v in &layer {
            new_layer.push(*v * values[i]);
        }
        layer = new_layer;
    }

    layer
}
