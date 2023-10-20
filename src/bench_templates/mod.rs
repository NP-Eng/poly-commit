use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{rand::Rng, test_rng};

/// type alias for DenseMultilinearExtension
pub type MLE<F> = DenseMultilinearExtension<F>;

use core::time::Duration;
use std::time::Instant;

use crate::{challenge::ChallengeGenerator, LabeledPolynomial, PolynomialCommitment};

use criterion::{BenchmarkId, Criterion};

/// Measure the time cost of {commit/open/verify} across a range of num_vars
pub fn bench_pcs_method<
    F: PrimeField,
    PCS: PolynomialCommitment<F, DenseMultilinearExtension<F>, PoseidonSponge<F>>,
>(
    c: &mut Criterion,
    range: Vec<usize>,
    msg: &str,
    method: impl Fn(&PCS::UniversalParams, usize) -> Duration,
) {
    let mut group = c.benchmark_group(msg);
    let rng = &mut test_rng();

    // Add for logarithmic scale (should yield linear plots)
    // let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    // group.plot_config(plot_config);

    for num_vars in range {
        // TODO if this takes too long and key trimming works, we might want to pull this out from the loop
        let pp = PCS::setup(1, Some(num_vars), rng).unwrap();

        group.bench_with_input(
            BenchmarkId::from_parameter(num_vars),
            &num_vars,
            |b, num_vars| {
                b.iter(|| method(&pp, *num_vars));
            },
        );
    }

    group.finish();
}

/// Report the time cost of a commitment
pub fn commit<
    F: PrimeField,
    PCS: PolynomialCommitment<F, DenseMultilinearExtension<F>, PoseidonSponge<F>>,
>(
    pp: &PCS::UniversalParams,
    num_vars: usize,
) -> Duration {
    // TODO create or pass? depends on the cost
    let rng = &mut test_rng();

    let (ck, _) = PCS::trim(&pp, 1, 1, None).unwrap();

    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None);

    let start = Instant::now();
    let (_, _) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    start.elapsed()
}

/// Report the size of a commitment
pub fn commitment_size<
    F: PrimeField,
    PCS: PolynomialCommitment<F, DenseMultilinearExtension<F>, PoseidonSponge<F>>,
>(
    num_vars: usize,
) -> usize {
    let rng = &mut test_rng();
    let pp = PCS::setup(1, Some(num_vars), rng).unwrap();

    let (ck, _) = PCS::trim(&pp, 1, 1, None).unwrap();

    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None);

    let (coms, _) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();

    coms[0].commitment().serialized_size(Compress::No)
}

/// Report the time cost of an opening
pub fn open<
    F: PrimeField,
    PCS: PolynomialCommitment<F, DenseMultilinearExtension<F>, PoseidonSponge<F>>,
>(
    pp: &PCS::UniversalParams,
    num_vars: usize,
) -> Duration {
    let rng = &mut test_rng();
    let (ck, _) = PCS::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None);

    let (coms, randomness) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = rand_mv_point(num_vars, rng);

    let start = Instant::now();
    let _ = PCS::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut ChallengeGenerator::new_univariate(&mut test_sponge()),
        &randomness,
        Some(rng),
    )
    .unwrap();
    start.elapsed()
}

/// Report the size of a proof
pub fn proof_size<
    F: PrimeField,
    PCS: PolynomialCommitment<F, DenseMultilinearExtension<F>, PoseidonSponge<F>>,
>(
    num_vars: usize,
) -> usize {
    let rng = &mut test_rng();
    let pp = PCS::setup(1, Some(num_vars), rng).unwrap();

    let (ck, _) = PCS::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None);

    let (coms, randomness) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = rand_mv_point(num_vars, rng);

    let proofs = PCS::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut ChallengeGenerator::new_univariate(&mut test_sponge()),
        &randomness,
        Some(rng),
    )
    .unwrap();

    let bproof: PCS::BatchProof = vec![proofs].into();

    bproof.serialized_size(Compress::No)
}

/// Report the time cost of a verification
pub fn verify<
    F: PrimeField,
    PCS: PolynomialCommitment<F, DenseMultilinearExtension<F>, PoseidonSponge<F>>,
>(
    pp: &PCS::UniversalParams,
    num_vars: usize,
) -> Duration {
    let rng = &mut test_rng();
    let (ck, vk) = PCS::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_ml_poly(num_vars, rng), None, None);

    let (coms, randomness) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = rand_mv_point(num_vars, rng);
    let claimed_eval = labeled_poly.evaluate(&point);
    let proof = PCS::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut ChallengeGenerator::new_univariate(&mut test_sponge()),
        &randomness,
        Some(rng),
    )
    .unwrap();

    let start = Instant::now();
    PCS::check(
        &vk,
        &coms,
        &point,
        [claimed_eval],
        &proof,
        &mut ChallengeGenerator::new_univariate(&mut test_sponge()),
        None,
    )
    .unwrap();
    start.elapsed()
}

/*************** Auxiliary functions ***************/

fn rand_ml_poly<F: PrimeField>(num_vars: usize, rng: &mut impl Rng) -> MLE<F> {
    MLE::rand(num_vars, rng)
}

fn rand_mv_point<F: PrimeField>(num_vars: usize, rng: &mut impl Rng) -> Vec<F> {
    (0..num_vars).map(|_| F::rand(rng)).collect()
}

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
