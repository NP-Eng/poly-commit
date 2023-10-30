use ark_ec::AffineRepr;
use ark_pcs_bench_templates::*;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_bn254::{Fr, G1Affine};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial as DenseUnivariatePoly;
use ark_poly_commit::ipa_pc::InnerProductArgPC;

use rand_chacha::ChaCha20Rng;

type UniPoly = DenseUnivariatePoly<Fr>;
type Sponge = PoseidonSponge<<EdwardsAffine as AffineRepr>::ScalarField>;

// Hyrax PCS over BN254
type Hyrax254 = HyraxPC<G1Affine, DenseMultilinearExtension<Fr>>;

fn rand_poly_hyrax<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> DenseMultilinearExtension<F> {
    DenseMultilinearExtension::rand(num_vars, rng)
}

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 20;

bench!(Hyrax254, rand_poly_hyrax);