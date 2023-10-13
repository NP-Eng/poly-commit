use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ec::AffineRepr;
use ark_ed_on_bls12_381::EdwardsAffine;
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::test_rng;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use crate::challenge::ChallengeGenerator;
use crate::hyrax::HyraxPC;

use crate::linear_codes::utils::test_sponge;
use crate::{LabeledPolynomial, PolynomialCommitment};

use crate::tests::*;

// The test structure is largely taken from the multilinear_ligero module
// inside this crate

// ****************** types ******************

type Fr = <EdwardsAffine as AffineRepr>::ScalarField;
type Hyrax = HyraxPC<EdwardsAffine, DenseMultilinearExtension<Fr>>;

// ******** auxiliary test functions ********

fn rand_poly<Fr: PrimeField>(
    _: usize, // degree
    num_vars: Option<usize>,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<Fr> {
    match num_vars {
        Some(n) => DenseMultilinearExtension::rand(n, rng),
        None => panic!("Must specify number of variables"),
    }
}

fn constant_poly<Fr: PrimeField>(
    _: usize, // degree
    num_vars: Option<usize>,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<Fr> {
    match num_vars {
        Some(0) => DenseMultilinearExtension::rand(0, rng),
        _ => panic!("Must specify number of variables: 0"),
    }
}

fn rand_point<F: PrimeField>(num_vars: Option<usize>, rng: &mut ChaCha20Rng) -> Vec<F> {
    match num_vars {
        Some(n) => (0..n).map(|_| F::rand(rng)).collect(),
        None => panic!("Must specify number of variables"),
    }
}

// ****************** tests ******************

#[test]
fn test_hyrax_construction() {
    // Desired number of variables (must be even!)
    let n = 8;

    let chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let pp = Hyrax::setup(1, Some(n), chacha).unwrap();

    let (ck, vk) = Hyrax::trim(&pp, 1, 1, None).unwrap();

    let l_poly = LabeledPolynomial::new(
        "test_poly".to_string(),
        rand_poly::<Fr>(0, Some(n), chacha),
        None,
        None,
    );

    let (c, rands) = Hyrax::commit(&ck, &[l_poly.clone()], Some(chacha)).unwrap();

    let point: Vec<Fr> = rand_point(Some(n), chacha);
    let value = l_poly.evaluate(&point);

    // Dummy argument
    let mut test_sponge = test_sponge::<Fr>();
    let mut challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
        ChallengeGenerator::new_univariate(&mut test_sponge);

    let proof = Hyrax::open(
        &ck,
        &[l_poly],
        &c,
        &point,
        &mut (challenge_generator.clone()),
        &rands,
        Some(chacha),
    )
    .unwrap();

    assert!(Hyrax::check(
        &vk,
        &c,
        &point,
        [value],
        &proof,
        &mut challenge_generator,
        Some(chacha),
    )
    .unwrap());
}

#[test]
fn hyrax_single_poly_test() {
    single_poly_test::<_, _, Hyrax, _>(
        Some(10),
        rand_poly::<Fr>,
        rand_point::<Fr>,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
    // single_poly_test::<_, _, HyraxPC<_, _>, _>(
    //     Some(10),
    //     rand_poly::<Fr381>,
    //     rand_point::<Fr381>,
    //     poseidon_sponge_for_test,
    // )
    // .expect("test failed for bls12-381");
}

#[test]
fn hyrax_constant_poly_test() {
    // use crate::tests::*;
    // single_poly_test::<_, _, LigeroPCS, _>(
    //     Some(10),
    //     constant_poly::<Fr>,
    //     rand_point::<Fr>,
    //     poseidon_sponge_for_test,
    // )
    // .expect("test failed for bls12-377");
    single_poly_test::<_, _, Hyrax, _>(
        Some(0),
        constant_poly::<Fr>,
        rand_point::<Fr>,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_full_end_to_end_test() {
    // full_end_to_end_test::<_, _, LigeroPCS, _>(
    //     Some(8),
    //     rand_poly::<Fr>,
    //     rand_point::<Fr>,
    //     poseidon_sponge_for_test,
    // )
    // .expect("test failed for bls12-377");
    full_end_to_end_test::<_, _, Hyrax, _>(
        Some(10),
        rand_poly::<Fr>,
        rand_point::<Fr>,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_single_equation_test() {
    // single_equation_test::<_, _, LigeroPCS, _>(
    //     Some(10),
    //     rand_poly::<Fr>,
    //     rand_point::<Fr>,
    //     poseidon_sponge_for_test,
    // )
    // .expect("test failed for bls12-377");
    single_equation_test::<_, _, Hyrax, _>(
        Some(6),
        rand_poly::<Fr>,
        rand_point::<Fr>,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_two_equation_test() {
    // two_equation_test::<_, _, LigeroPCS, _>(
    //     Some(5),
    //     rand_poly::<Fr>,
    //     rand_point::<Fr>,
    //     poseidon_sponge_for_test,
    // )
    // .expect("test failed for bls12-377");
    two_equation_test::<_, _, Hyrax, _>(
        Some(10),
        rand_poly::<Fr>,
        rand_point::<Fr>,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_full_end_to_end_equation_test() {
    // full_end_to_end_equation_test::<_, _, LigeroPCS, _>(
    //     Some(5),
    //     rand_poly::<Fr>,
    //     rand_point::<Fr>,
    //     poseidon_sponge_for_test,
    // )
    println!("Finished bls12-377");
    full_end_to_end_equation_test::<_, _, Hyrax, _>(
        Some(8),
        rand_poly::<Fr>,
        rand_point::<Fr>,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}
