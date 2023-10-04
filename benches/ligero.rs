use std::time::SystemTime;

use ark_bls12_377::Fq;
use ark_bls12_377::Fr;
use ark_bls12_381::Fr as Fr381;
use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_crypto_primitives::{
    crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
    sponge::poseidon::PoseidonSponge,
};
use ark_ff::{Field, PrimeField};
use ark_poly::{
    domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
    EvaluationDomain, Polynomial,
};
use ark_poly_commit::LabeledCommitment;
use ark_poly_commit::{
    challenge::ChallengeGenerator,
    ligero::{Ligero, LigeroPCUniversalParams},
    LabeledPolynomial, PolynomialCommitment,
};
use ark_std::rand::Rng;
use ark_std::test_rng;
use ark_std::UniformRand;
use blake2::Blake2s256;
use criterion::{criterion_group, criterion_main, Criterion};
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

struct MerkleTreeParams;
type LeafH = Sha256;
type CompressH = Sha256;

impl Config for MerkleTreeParams {
    type Leaf = [u8];

    type LeafDigest = <LeafH as CRHScheme>::Output;
    type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

type MTConfig = MerkleTreeParams;
type UniPoly = DensePolynomial<Fr>;
type Sponge = PoseidonSponge<Fr>;
type PC<F, C, D, S, P> = Ligero<F, C, D, S, P>;
type LigeroPCS = PC<Fr, MTConfig, Blake2s256, Sponge, UniPoly>;
// TODO are we going to bench with various fields?
type _LigeroPcsF<F> = PC<F, MTConfig, Blake2s256, Sponge, DensePolynomial<F>>;

fn rand_poly<F: PrimeField>(
    degree: usize,
    _: Option<usize>,
    rng: &mut impl Rng,
) -> DensePolynomial<F> {
    DensePolynomial::rand(degree, rng)
}

const SAMPLES: usize = 100;

fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut v = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();

        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}

fn commit(c: &mut Criterion) {
    // degree is 18 like in Jellyfish Multilinear KZG
    let degree = 18;

    let rng = &mut test_rng();
    let pp = LigeroPCS::setup(degree, None, rng).unwrap();
    let (ck, _) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| {
            LabeledPolynomial::new("test".to_string(), rand_poly(degree, None, rng), None, None)
        })
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labaled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    c.bench_function("Ligero Commit", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            let (_, _) = LigeroPCS::commit(&ck, [labaled_poly_refs[i]], None).unwrap();
        })
    });
}

fn open(c: &mut Criterion) {
    // degree is 18 like in Jellyfish Multilinear KZG
    let degree = 18;

    let rng = &mut test_rng();
    let pp = LigeroPCS::setup(degree, None, rng).unwrap();
    let (ck, _) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| {
            LabeledPolynomial::new("test".to_string(), rand_poly(degree, None, rng), None, None)
        })
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labaled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    let commitments: Vec<(Vec<_>, _)> = (0..SAMPLES)
        .map(|i| {
            let (c, r) = LigeroPCS::commit(&ck, [labaled_poly_refs[i]], None).unwrap();
            (c, r)
        })
        .collect();
    let challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
        ChallengeGenerator::new_univariate(&mut test_sponge());

    let points: Vec<_> = (0..SAMPLES).map(|_| Fr::rand(&mut test_rng())).collect();

    c.bench_function("Ligero Open", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            LigeroPCS::open(
                &ck,
                [labaled_poly_refs[i]],
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

criterion_group! {
    name = ligero_benches;
    config = Criterion::default();
    targets =
        commit,
        open,
}

criterion_main!(ligero_benches);
