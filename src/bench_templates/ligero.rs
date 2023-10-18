use core::num;

use ark_bls12_377::Fr;
use ark_crypto_primitives::{
    crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
    sponge::{
        poseidon::{PoseidonConfig, PoseidonSponge},
        CryptographicSponge,
    },
};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_poly_commit::{
    challenge::ChallengeGenerator,
    linear_codes::{FieldToBytesColHasher, LeafIdentityHasher, LinearCodePCS, MultilinearLigero},
    LabeledPolynomial, PolynomialCommitment,
    hyrax::HyraxPC
};
use ark_std::rand::Rng;
use ark_std::test_rng;
use ark_std::UniformRand;
use blake2::Blake2s256;
use criterion::{criterion_group, criterion_main, Criterion};
struct MerkleTreeParams;
type LeafH = LeafIdentityHasher;
type CompressH = Sha256;

// ******** Types

impl Config for MerkleTreeParams {
    type Leaf = Vec<u8>;

    type LeafDigest = <LeafH as CRHScheme>::Output;
    type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

type MTConfig = MerkleTreeParams;
type Sponge = PoseidonSponge<Fr>;
type ColHasher = FieldToBytesColHasher<Fr, Blake2s256>;
type LigeroPCS = LinearCodePCS<
    MultilinearLigero<Fr, MTConfig, Sponge, MLE, ColHasher>,
    Fr,
    MLE,
    Sponge,
    MTConfig,
    ColHasher,
>;

fn commit_ligero(c: &mut Criterion) {
    let num_vars = 18;

    let rng = &mut test_rng();
    let pp = LigeroPCS::setup(num_vars, None, rng).unwrap();
    let (ck, _) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None))
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labeled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    c.bench_function("Ligero Commit", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            let (_, _) = LigeroPCS::commit(&ck, [labeled_poly_refs[i]], None).unwrap();
        })
    });
}

fn open_ligero(c: &mut Criterion) {
    let num_vars = 18;

    let rng = &mut test_rng();
    let pp = LigeroPCS::setup(num_vars, None, rng).unwrap();
    let (ck, _) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None))
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labeled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    let commitments: Vec<(Vec<_>, _)> = (0..SAMPLES)
        .map(|i| {
            let (c, r) = LigeroPCS::commit(&ck, [labeled_poly_refs[i]], None).unwrap();
            (c, r)
        })
        .collect();
    let challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
        ChallengeGenerator::new_univariate(&mut test_sponge());

    let points: Vec<_> = (0..SAMPLES)
        .map(|_| rand_point(Some(num_vars), &mut test_rng()))
        .collect();

    c.bench_function("Ligero Open", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            let _commitment = LigeroPCS::open(
                &ck,
                [labeled_poly_refs[i]],
                &commitments[i].0,
                &points[i],
                &mut (challenge_generator.clone()),
                &commitments[i].1,
                None,
            )
            .unwrap();
        })
    });
}

fn verify_ligero(c: &mut Criterion) {
    let num_vars = 18;

    let rng = &mut test_rng();
    let pp = LigeroPCS::setup(num_vars, None, rng).unwrap();
    let (ck, vk) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None))
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labeled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    let commitments: Vec<(Vec<_>, _)> = (0..SAMPLES)
        .map(|i| {
            let (c, r) = LigeroPCS::commit(&ck, [labeled_poly_refs[i]], None).unwrap();
            (c, r)
        })
        .collect();
    let challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
        ChallengeGenerator::new_univariate(&mut test_sponge());

    let points: Vec<_> = (0..SAMPLES)
        .map(|_| rand_point(Some(num_vars), &mut test_rng()))
        .collect();

    let claimed_evals = (0..SAMPLES)
        .map(|i| (labeled_polys[i].evaluate(&points[i])))
        .collect::<Vec<_>>();

    // instead of benching open, we create a list of proofs
    let proofs: Vec<_> = (0..SAMPLES)
        .map(|i| {
            LigeroPCS::open(
                &ck,
                [labeled_poly_refs[i]],
                &commitments[i].0,
                &points[i],
                &mut (challenge_generator.clone()),
                &commitments[i].1,
                None,
            )
            .unwrap()
        })
        .collect();

    c.bench_function("Ligero Verify", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;

            // we check the proof
            LigeroPCS::check(
                &vk,
                &commitments[i].0,
                &points[i],
                [claimed_evals[i]],
                &proofs[i],
                &mut (challenge_generator.clone()),
                None,
            ).unwrap();
        })
    });
}
