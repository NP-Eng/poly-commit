use core::time::Duration;
use std::time::Instant;

use crate::{
    challenge::ChallengeGenerator, hyrax::HyraxPC, LabeledPolynomial, PolynomialCommitment,
};
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ec::AffineRepr;
use ark_poly::DenseMultilinearExtension;
use ark_std::test_rng;
use criterion::Criterion;

use super::{rand_ml_poly, rand_mv_point, test_sponge, SAMPLES};

/************ Hyrax types ************/
//type Hyrax377 = HyraxPC<G1Affine, DenseMultilinearExtension<Fr>>;
type Hyrax<G> = HyraxPC<G, DenseMultilinearExtension<<G as AffineRepr>::ScalarField>>;

/************ Hyrax functions ************/

/// Measure the time cost of Hyrax commitments
pub fn commit_hyrax<G: AffineRepr>(c: &mut Criterion, num_vars: usize, msg: &str) {
    let rng = &mut test_rng();
    let pp = Hyrax::<G>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Hyrax::<G>::trim(&pp, 1, 1, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| {
            LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None)
        })
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labeled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    c.bench_function(msg, |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            let (_, _) = Hyrax::<G>::commit(&ck, [labeled_poly_refs[i]], Some(rng)).unwrap();
        })
    });
}

/// Measure the time cost of Hyrax commitments over a range of numbers of variables
pub fn commit_hyrax_custom<G: AffineRepr>(num_vars: usize, msg: &str) -> Duration {
    let rng = &mut test_rng();
    let pp = Hyrax::<G>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Hyrax::<G>::trim(&pp, 1, 1, None).unwrap();

    let labeled_poly = LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None);

    let start = Instant::now();
    let (_, _) = Hyrax::<G>::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    start.elapsed()
}

/// Measure the time cost of Hyrax openings
pub fn open_hyrax<G: AffineRepr>(c: &mut Criterion, num_vars: usize, msg: &str) {
    let rng = &mut test_rng();
    let pp = Hyrax::<G>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Hyrax::<G>::trim(&pp, 1, 1, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| {
            LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None)
        })
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labeled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    let commitments: Vec<(Vec<_>, _)> = (0..SAMPLES)
        .map(|i| {
            let (c, r) = Hyrax::<G>::commit(&ck, [labeled_poly_refs[i]], Some(rng)).unwrap();
            (c, r)
        })
        .collect();

    let challenge_generator: ChallengeGenerator<G::ScalarField, PoseidonSponge<G::ScalarField>> =
        ChallengeGenerator::new_univariate(&mut test_sponge());

    let points: Vec<_> = (0..SAMPLES)
        .map(|_| rand_mv_point(num_vars, &mut test_rng()))
        .collect();

    c.bench_function(msg, |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            let _commitment = Hyrax::<G>::open(
                &ck,
                [labeled_poly_refs[i]],
                &commitments[i].0,
                &points[i],
                &mut (challenge_generator.clone()),
                &commitments[i].1,
                Some(rng),
            )
            .unwrap();
        })
    });
}

/// Measure the time cost of Hyrax proof verification
pub fn verify_hyrax<G: AffineRepr>(c: &mut Criterion, num_vars: usize, msg: &str) {
    let rng = &mut test_rng();
    let pp = Hyrax::<G>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, vk) = Hyrax::<G>::trim(&pp, 1, 1, None).unwrap();

    let labeled_polys = (0..SAMPLES)
        .map(|_| {
            LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None)
        })
        .collect::<Vec<_>>();

    // this is a little ugly, but ideally we want to avoid cloning inside the benchmark. Therefore we keep `labeled_polys` in scope, and just commit to references to it.
    let labeled_poly_refs = labeled_polys.iter().map(|p| p).collect::<Vec<_>>();

    let commitments: Vec<(Vec<_>, _)> = (0..SAMPLES)
        .map(|i| {
            let (c, r) = Hyrax::<G>::commit(&ck, [labeled_poly_refs[i]], Some(rng)).unwrap();
            (c, r)
        })
        .collect();

    let challenge_generator: ChallengeGenerator<G::ScalarField, PoseidonSponge<G::ScalarField>> =
        ChallengeGenerator::new_univariate(&mut test_sponge());

    let points: Vec<_> = (0..SAMPLES)
        .map(|_| rand_mv_point(num_vars, &mut test_rng()))
        .collect();

    let claimed_evals = (0..SAMPLES)
        .map(|i| (labeled_polys[i].evaluate(&points[i])))
        .collect::<Vec<_>>();

    let proofs: Vec<_> = (0..SAMPLES)
        .map(|i| {
            Hyrax::<G>::open(
                &ck,
                [labeled_poly_refs[i]],
                &commitments[i].0,
                &points[i],
                &mut (challenge_generator.clone()),
                &commitments[i].1,
                Some(rng),
            )
            .unwrap()
        })
        .collect();

    c.bench_function(msg, |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;

            // we check the proof
            Hyrax::<G>::check(
                &vk,
                &commitments[i].0,
                &points[i],
                [claimed_evals[i]],
                &proofs[i],
                &mut (challenge_generator.clone()),
                None,
            )
            .unwrap();
        })
    });
}
