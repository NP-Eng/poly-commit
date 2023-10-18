#![cfg(feature = "benches")]
use core::num;

use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use criterion::{criterion_group, criterion_main, Criterion};

use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_std::{rand::Rng, test_rng};
use criterion::BenchmarkId;

use ark_poly_commit::bench_templates::hyrax::{commit_hyrax, commit_hyrax_range};

use ark_bls12_377::G1Affine;
use ark_ed_on_bls12_381::EdwardsAffine;

/*************** Instantiating target functions ***************/
fn commit_hyrax_n12_c377(c: &mut Criterion) {
    commit_hyrax::<G1Affine>(c, 12, "commit_hyrax_n12_c377");
}
fn commit_hyrax_n14_c377(c: &mut Criterion) {
    commit_hyrax::<G1Affine>(c, 14, "commit_hyrax_n14_c377");
}

fn commit_hyrax_n12_c381(c: &mut Criterion) {
    commit_hyrax::<EdwardsAffine>(c, 12, "commit_hyrax_n12_c381");
}
fn commit_hyrax_n14_c381(c: &mut Criterion) {
    commit_hyrax::<EdwardsAffine>(c, 14, "commit_hyrax_n14_c381");
}

fn commit_hyrax_range_c377(c: &mut Criterion) {
    commit_hyrax_range::<G1Affine>(c, (10..20).step_by(2).collect(), "commit_hyrax_range_c377");
}

fn commit_hyrax_range_c381(c: &mut Criterion) {
    commit_hyrax_range::<EdwardsAffine>(
        c,
        (10..20).step_by(2).collect(),
        "commit_hyrax_range_c381",
    );
}

criterion_group!(hyrax_benches_range_377, commit_hyrax_range_c377);
criterion_group!(hyrax_benches_range_381, commit_hyrax_range_c381);

/*************** Benchmarks ***************/

// criterion_group! {
//     name = mixed_benches;
//     config = Criterion::default();
//     targets =
//         commit_ligero,
//         commit_hyrax,
//         open_ligero,
//         open_hyrax,
//         verify_ligero,
//         verify_hyrax,
// }

criterion_group! {
    name = hyrax_benches;
    config = Criterion::default();
    targets =
        commit_hyrax_n12_c377,
        commit_hyrax_n14_c377,
        commit_hyrax_n12_c381,
        commit_hyrax_n14_c381,
}

criterion_main!(hyrax_benches_range_377);
