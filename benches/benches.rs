#![cfg(feature = "benches")]
use ark_ec::AffineRepr;
use ark_poly::DenseMultilinearExtension;
use blake2::Blake2s256;
use criterion::{criterion_group, criterion_main, Criterion};

use ark_crypto_primitives::{
    crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
    sponge::poseidon::PoseidonSponge,
};

use ark_poly_commit::{
    bench_templates::{bench_pcs_method, commit, open, verify, MLE},
    hyrax::HyraxPC,
    linear_codes::{FieldToBytesColHasher, LeafIdentityHasher, LinearCodePCS, MultilinearLigero, MultilinearBrakedown},
};

use ark_bls12_381::{Fr as Fr381, G1Affine as G1Affine381};
use ark_bn254::{Fr as Fr254, G1Affine as G1Affine254};

// Hyrax type alias
type Hyrax<G> = HyraxPC<G, DenseMultilinearExtension<<G as AffineRepr>::ScalarField>>;

struct MerkleTreeParams;
type LeafH = LeafIdentityHasher;
type CompressH = Sha256;
impl Config for MerkleTreeParams {
    type Leaf = Vec<u8>;

    type LeafDigest = <LeafH as CRHScheme>::Output;
    type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

type MTConfig = MerkleTreeParams;
type Sponge<F> = PoseidonSponge<F>;
type ColHasher<F> = FieldToBytesColHasher<F, Blake2s256>;

// Ligero type alias
type Ligero<F> = LinearCodePCS<
    MultilinearLigero<F, MTConfig, Sponge<F>, MLE<F>, ColHasher<F>>,
    F,
    MLE<F>,
    Sponge<F>,
    MTConfig,
    ColHasher<F>,
>;

// Brakedown type alias
type Brakedown<F> = LinearCodePCS<
    MultilinearBrakedown<F, MTConfig, Sponge<F>, MLE<F>, ColHasher<F>>,
    F,
    MLE<F>,
    Sponge<F>,
    MTConfig,
    ColHasher<F>,
>;

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 20;

/*************** Instantiating target functions ***************/
fn hyrax_bls12_381(c: &mut Criterion) {
    bench_pcs_method::<_, Hyrax<G1Affine381>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "commit_hyrax_range_BLS12_381",
        commit::<_, Hyrax<G1Affine381>>,
    );
    bench_pcs_method::<_, Hyrax<G1Affine381>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "open_hyrax_range_BLS12_381",
        open::<_, Hyrax<G1Affine381>>,
    );

    bench_pcs_method::<_, Hyrax<G1Affine381>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "verify_hyrax_range_BLS12_381",
        verify::<_, Hyrax<G1Affine381>>,
    );
}

fn hyrax_bn254(c: &mut Criterion) {
    bench_pcs_method::<_, Hyrax<G1Affine254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "commit_hyrax_range_BN_254",
        commit::<_, Hyrax<G1Affine254>>,
    );
    bench_pcs_method::<_, Hyrax<G1Affine254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "open_hyrax_range_BN_254",
        open::<_, Hyrax<G1Affine254>>,
    );

    bench_pcs_method::<_, Hyrax<G1Affine254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "verify_hyrax_range_BN_254",
        verify::<_, Hyrax<G1Affine254>>,
    );
}

fn ligero_bls12_381(c: &mut Criterion) {
    bench_pcs_method::<_, Ligero<Fr381>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "commit_ligero_range_BLS12_381",
        commit::<_, Ligero<Fr381>>,
    );
    bench_pcs_method::<_, Ligero<Fr381>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "open_ligero_range_BLS12_381",
        open::<_, Ligero<Fr381>>,
    );

    bench_pcs_method::<_, Ligero<Fr381>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "verify_ligero_range_BLS12_381",
        verify::<_, Ligero<Fr381>>,
    );
}

fn ligero_bn254(c: &mut Criterion) {
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

    bench_pcs_method::<_, Ligero<Fr254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "verify_ligero_range_BN_254",
        verify::<_, Ligero<Fr254>>,
    );
}

fn brakedown_bn254(c: &mut Criterion) {
    bench_pcs_method::<_, Brakedown<Fr254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "commit_brakedown_range_BN_254",
        commit::<_, Brakedown<Fr254>>,
    );
    bench_pcs_method::<_, Brakedown<Fr254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "open_brakedown_range_BN_254",
        open::<_, Brakedown<Fr254>>,
    );

    bench_pcs_method::<_, Brakedown<Fr254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
        "verify_brakedown_range_BN_254",
        verify::<_, Brakedown<Fr254>>,
    );
}

criterion_group! {
    name = hyrax_benches;
    config = Criterion::default();
    targets =
        hyrax_bls12_381,
        hyrax_bn254
}

criterion_group! {
    name = ligero_benches;
    config = Criterion::default();
    targets =
        ligero_bls12_381,
        ligero_bn254
}

criterion_group! {
    name = brakedown_benches;
    config = Criterion::default();
    targets =
        brakedown_bn254,
}

criterion_main!(hyrax_benches, ligero_benches, brakedown_benches);
