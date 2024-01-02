use ark_pcs_bench_templates::*;
use ark_poly::DenseMultilinearExtension;

use ark_bn254::{Fr, G1Affine};
use ark_ff::{BigInteger, PrimeField};
use ark_poly_commit::hyrax::HyraxPC;

use rand_chacha::ChaCha20Rng;

// Hyrax PCS over BN254
type Hyrax254 = HyraxPC<G1Affine, DenseMultilinearExtension<Fr>>;

fn rand_poly_hyrax<F: PrimeField>(
    num_vars: usize,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<F> {
    let max_bits: usize = 30;
    let num_bits = F::MODULUS_BIT_SIZE as usize;
    let small_scalars = (0..(1 << num_vars))
        .map(|_| {
            let s = F::rand(rng).into_bigint();
            let mut bits = s.to_bits_le();
            bits.truncate(max_bits);
            bits.resize(num_bits, false);
            let bigint = F::BigInt::from_bits_le(&bits);
            F::from_bigint(bigint).unwrap()
        })
        .collect::<Vec<_>>();

    DenseMultilinearExtension::from_evaluations_vec(num_vars, small_scalars)
}

fn rand_point_hyrax<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
    (0..num_vars).map(|_| F::rand(rng)).collect()
}

const MIN_NUM_VARS: usize = 32;
const MAX_NUM_VARS: usize = 34;

bench!(Hyrax254, rand_poly_hyrax, rand_point_hyrax);
