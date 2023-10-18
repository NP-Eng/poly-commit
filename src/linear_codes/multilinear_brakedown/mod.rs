use super::{BrakedownPCParams, LinearEncode};
use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{merkle_tree::Config, sponge::CryptographicSponge};
use ark_ff::{Field, PrimeField};
use ark_poly::{MultilinearExtension, Polynomial};
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;

mod tests;

/// The univariate Brakedown polynomial commitment scheme based on [[Brakedown]][bd].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [bd]: https://eprint.iacr.org/2021/1043.pdf
pub struct MultilinearBrakedown<
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    H: CRHScheme,
> {
    _phantom: PhantomData<(F, C, S, P, H)>,
}

impl<F, C, S, P, H> LinearEncode<F, C, P, H> for MultilinearBrakedown<F, C, S, P, H>
where
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    <P as Polynomial<F>>::Point: Into<Vec<F>>,
    H: CRHScheme,
{
    type LinCodePCParams = BrakedownPCParams<F, C, H>;

    fn setup<R: RngCore>(
        _max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
        leaf_hash_params: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
        col_hash_params: H::Parameters,
    ) -> Self::LinCodePCParams {
        Self::LinCodePCParams::new(
            rng,
            128,
            (1, 5),
            (41, 500),
            (41, 25),
            30,
            1 << num_vars.unwrap(),
            true,
            leaf_hash_params,
            two_to_one_params,
            col_hash_params,
        )
    }

    fn encode(msg: &[F], pp: &Self::LinCodePCParams) -> Vec<F> {
        assert!(msg.len() == pp.n); // TODO Make it error
        let cw_len = pp.codeword_len();
        let mut cw = vec![F::zero(); cw_len];
        cw[..msg.len()].copy_from_slice(msg);

        // istart and iend are two arrays which keep the indices of the encoding in each step
        let start = pp
            .a_dims
            .iter()
            .scan(0, |acc, &(x, _, _)| {
                *acc += x;
                Some(*acc)
            })
            .collect::<Vec<_>>();
        let end = pp
            .b_dims
            .iter()
            .scan(cw_len, |acc, &(_, x, _)| {
                *acc -= x;
                Some(*acc)
            })
            .collect::<Vec<_>>();

        // Multiply by matrices A
        for (i, &s) in start.iter().enumerate() {
            let src = &pp.a_mats[i].row_mul(&cw[s - pp.a_dims[i].0..s]);
            cw[s..s + pp.a_dims[i].1].copy_from_slice(src);
        }

        // RS encode the last one
        let rss = if let Some(s) = start.last() { *s } else { 0 };
        let rsie = if let Some(e) = pp.a_dims.last() {
            rss + e.1
        } else {
            pp.n
        };
        let rsoe = if let Some(e) = end.last() { *e } else { cw_len };
        naive_reed_solomon(&mut cw, rss, rsie, rsoe);

        // Come back
        for (i, (&s, e)) in start.iter().zip(end).enumerate() {
            let src = &pp.b_mats[i].row_mul(&cw[s..e]);
            cw[e..e + pp.b_dims[i].1].copy_from_slice(src);
        }
        cw.to_vec()
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
}

// This RS encoding is on points 1, ..., oe - s without relying on FFTs
fn naive_reed_solomon<F: Field>(cw: &mut [F], s: usize, ie: usize, oe: usize) {
    let mut res = vec![F::zero(); oe - s];
    let mut x = F::one();
    for r in res.iter_mut() {
        for j in (s..ie).rev() {
            *r *= x;
            *r += cw[j];
        }
        x += F::one();
    }
    cw[s..oe].copy_from_slice(&res);
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
