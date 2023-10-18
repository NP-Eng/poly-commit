#![cfg(feature = "benches")]
use std::{borrow::Borrow, marker::PhantomData};

use ark_ff::PrimeField;
use ark_std::{rand::RngCore, test_rng};
use criterion::{criterion_group, criterion_main, Criterion};

use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{Config, IdentityDigestConverter},
    sponge::{
        poseidon::{PoseidonConfig, PoseidonSponge},
        Absorb, CryptographicSponge,
    },
    Error,
};

use ark_poly_commit::{
    bench_templates::{bench_pcs_method, commit, open, verify, MLE},
    linear_codes::{LinearCodePCS, MultilinearLigero},
};

/// Generate default parameters for alpha = 17, state-size = 8
///
/// WARNING: This poseidon parameter is not secure. Please generate
/// your own parameters according the field you use.
pub fn poseidon_parameters_for_test<F: PrimeField>() -> PoseidonConfig<F> {
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

use ark_bn254::Fr as Fr254;

// We introduce the wrapper only for the purpose of `setup` function not panicing with unimplemented
struct PoseidonWrapper<F>(PhantomData<F>);

impl<F: PrimeField + Absorb> CRHScheme for PoseidonWrapper<F> {
    type Input = [F];
    type Output = F;
    type Parameters = PoseidonConfig<F>;

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(poseidon_parameters_for_test())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        parameters: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, Error> {
        let input = input.borrow();

        let mut sponge = PoseidonSponge::new(parameters);
        sponge.absorb(&input);
        let res = sponge.squeeze_field_elements::<F>(1);
        Ok(res[0])
    }
}

impl<F: PrimeField + Absorb> TwoToOneCRHScheme for PoseidonWrapper<F> {
    type Input = F;
    type Output = F;
    type Parameters = PoseidonConfig<F>;

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(poseidon_parameters_for_test())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        parameters: &Self::Parameters,
        left_input: T,
        right_input: T,
    ) -> Result<Self::Output, Error> {
        let left_input = left_input.borrow();
        let right_input = right_input.borrow();

        let mut sponge = PoseidonSponge::new(parameters);
        sponge.absorb(left_input);
        sponge.absorb(right_input);
        let res = sponge.squeeze_field_elements::<F>(1);
        Ok(res[0])
    }

    fn compress<T: Borrow<Self::Output>>(
        parameters: &Self::Parameters,
        left_input: T,
        right_input: T,
    ) -> Result<Self::Output, Error> {
        let left_input = left_input.borrow();
        let right_input = right_input.borrow();

        let mut sponge = PoseidonSponge::new(parameters);
        sponge.absorb(left_input);
        sponge.absorb(right_input);
        let res = sponge.squeeze_field_elements::<F>(1);
        Ok(res[0])
    }
}

struct LeafIdentityHasher<F>(PhantomData<F>);

impl<F: PrimeField> CRHScheme for LeafIdentityHasher<F> {
    type Input = F;
    type Output = F;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, Error> {
        Ok(*input.borrow())
    }
}
struct MerkleTreeParams<F>(PhantomData<F>);

impl<F: PrimeField + Absorb> Config for MerkleTreeParams<F> {
    type Leaf = F;

    type LeafDigest = <LeafIdentityHasher<F> as CRHScheme>::Output;
    type LeafInnerDigestConverter = IdentityDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH<F> as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafIdentityHasher<F>;
    type TwoToOneHash = CompressH<F>;
}

type MTConfig<F> = MerkleTreeParams<F>;
type Sponge<F> = PoseidonSponge<F>;

type ColHasher<F> = PoseidonWrapper<F>;
type CompressH<F> = PoseidonWrapper<F>;
type Ligero<F> = LinearCodePCS<
    MultilinearLigero<F, MTConfig<F>, Sponge<F>, MLE<F>>,
    F,
    MLE<F>,
    Sponge<F>,
    MTConfig<F>,
    ColHasher<F>,
>;

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 20;

fn ligero_bn254(c: &mut Criterion) {
    // only commit and open; verify is done in-circuit!
    bench_pcs_method::<_, Ligero<Fr254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "commit_ligero_range_BN_254",
        commit::<_, Ligero<Fr254>>,
    );
    bench_pcs_method::<_, Ligero<Fr254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "open_ligero_range_BN_254",
        open::<_, Ligero<Fr254>>,
    );
}

criterion_group! {
    name = ligero_benches;
    config = Criterion::default();
    targets =
        ligero_bn254,
}

criterion_main!(ligero_benches);
