//#![cfg(feature = "benches")]

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{sponge::poseidon::PoseidonSponge, crh::sha256::Sha256};
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_poly_commit::linear_codes::{FieldToBytesColHasher, LinearCodePCS, MultilinearLigero};
use ark_poly_commit::{hyrax::HyraxPC, PolynomialCommitment, LabeledPolynomial, challenge::ChallengeGenerator, linear_codes::LeafIdentityHasher};
use ark_std::UniformRand;
use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, CryptographicSponge};
use ark_crypto_primitives::merkle_tree::{Config, ByteDigestConverter};

use ark_bls12_381::{Fr as Fr381, G1Affine as G1Affine381};
use ark_bn254::{Fr as Fr254, G1Affine as G1Affine254};
use ark_std::test_rng;
use blake2::Blake2s256;

use std::mem::size_of_val;

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 22;


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
type Ligero<F> = LinearCodePCS<
    MultilinearLigero<F, MTConfig, Sponge<F>, DenseMultilinearExtension<F>, ColHasher<F>>,
    F,
    DenseMultilinearExtension<F>,
    Sponge<F>,
    MTConfig,
    ColHasher<F>,
>;

/********************** Auxiliary functions *********************/

pub(crate) fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {

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

/********************** Size functions *********************/

fn hyrax_commitment_size<G: AffineRepr>(num_vars: usize) -> usize {

    let rng = &mut test_rng();
    let pp = Hyrax::<G>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Hyrax::<G>::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly = LabeledPolynomial::new("test".to_string(), 
        DenseMultilinearExtension::<G::ScalarField>::rand(num_vars, rng), None, None);

    let (coms, _) = Hyrax::<G>::commit(&ck, [&labeled_poly], Some(rng)).unwrap();

    size_of_val(&*coms[0].commitment().row_coms)
}

fn hyrax_proof_size<G: AffineRepr>(num_vars: usize) -> usize {

    let rng = &mut test_rng();
    let pp = Hyrax::<G>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Hyrax::<G>::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly = LabeledPolynomial::new("test".to_string(), 
        DenseMultilinearExtension::<G::ScalarField>::rand(num_vars, rng), None, None);

    let (coms, randomness) = Hyrax::<G>::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = (0..num_vars).map(|_| G::ScalarField::rand(rng)).collect();

    let proofs = Hyrax::<G>::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut ChallengeGenerator::new_univariate(&mut test_sponge()),
        &randomness,
        Some(rng),
    )
    .unwrap();

    let proof = proofs[0].clone();

    size_of_val(&proof) - size_of_val(&proof.z) + size_of_val(&*proof.z)
}

fn ligero_commitment_size<F: PrimeField>(num_vars: usize) -> usize {

    let rng = &mut test_rng();
    let pp = Ligero::<F>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Ligero::<F>::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly = LabeledPolynomial::new("test".to_string(), 
        DenseMultilinearExtension::<F>::rand(num_vars, rng), None, None);

    let (coms, _) = Ligero::<F>::commit(&ck, [&labeled_poly], Some(rng)).unwrap();

    size_of_val(&*coms[0].commitment().root)
}

fn ligero_proof_size<F: PrimeField>(num_vars: usize) -> usize {

    let rng = &mut test_rng();
    let pp = Ligero::<F>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, _) = Ligero::<F>::trim(&pp, 1, 1, None).unwrap();
    let labeled_poly = LabeledPolynomial::new("test".to_string(), 
        DenseMultilinearExtension::<F>::rand(num_vars, rng), None, None);

    let (coms, randomness) = Ligero::<F>::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = (0..num_vars).map(|_| F::rand(rng)).collect();

    let proofs = Ligero::<F>::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut ChallengeGenerator::new_univariate(&mut test_sponge()),
        &randomness,
        Some(rng),
    )
    .unwrap();

    let proof = proofs[0].clone();

    let mut size = 0;

    for p in proof.opening.paths {
        // TODO: probably incorrect
        size += size_of_val(&p);
    }

    for f in proof.opening.v {
        size += size_of_val(&f);
    }
    
    for col in proof.opening.columns {
        for f in col {
            size += size_of_val(&f);
        }
    }

    size
}

/********************** Main *********************/

fn main() {
    
    println!("\n---------------- Commitment size ----------------");

    println!("\nHyrax on BLS12-381: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, hyrax_commitment_size::<G1Affine381>(num_vars));
    }

    println!("\nHyrax on BN-254: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, hyrax_commitment_size::<G1Affine254>(num_vars));
    }

    println!("\nLigero on BLS12-381::Fr: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, ligero_commitment_size::<Fr381>(num_vars));
    }

    println!("\nLigero on BN-254::Fr: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, ligero_commitment_size::<Fr254>(num_vars));
    }

    println!("\n---------------- Proof size ----------------");

    println!("\nHyrax on BLS12-381: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, hyrax_proof_size::<G1Affine381>(num_vars));
    }

    println!("\nHyrax on BN-254: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, hyrax_proof_size::<G1Affine254>(num_vars));
    }

    println!("\nLigero on BLS12-381: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, ligero_proof_size::<Fr381>(num_vars));
    }

    println!("\nLigero on BN-254: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, ligero_proof_size::<Fr254>(num_vars));
    }
}
