use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{rand::Rng, test_rng};

// /// Auxiliary functions for benchmarking the Ligero PCS
// pub mod ligero;

/// Auxiliary functions for benchmarking the Hyrax PCS
pub mod hyrax;

type MLE<F> = DenseMultilinearExtension<F>;

pub(crate) const SAMPLES: usize = 100;

/*************** Auxiliary functions ***************/

fn rand_ml_poly<F: PrimeField>(num_vars: usize, rng: &mut impl Rng) -> MLE<F> {
    MLE::rand(num_vars, rng)
}

fn rand_uv_point<F: PrimeField>(rng: &mut impl Rng) -> F {
    F::rand(rng)
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
