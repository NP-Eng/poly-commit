#![cfg(feature = "benches")]

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{sponge::poseidon::PoseidonSponge, crh::sha256::Sha256};
use ark_ec::AffineRepr;
use ark_poly::DenseMultilinearExtension;
use ark_poly_commit::bench_templates::{commitment_size, proof_size};
use ark_poly_commit::linear_codes::{FieldToBytesColHasher, LinearCodePCS, MultilinearLigero};
use ark_poly_commit::{hyrax::HyraxPC, linear_codes::LeafIdentityHasher};
use ark_crypto_primitives::merkle_tree::{Config, ByteDigestConverter};

use ark_bls12_381::{Fr as Fr381, G1Affine as G1Affine381};
use ark_bn254::{Fr as Fr254, G1Affine as G1Affine254};
use blake2::Blake2s256;


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

/********************** Main *********************/

fn main() {
    
    println!("\n---------------- Commitment size ----------------");

    println!("\nHyrax on BLS12-381: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tser: {} B", commitment_size::<_, Hyrax<G1Affine381>>(num_vars));
    }

    println!("\nLigero on BLS12-381::Fr: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, commitment_size::<_, Ligero<Fr381>>(num_vars));
    }

    println!("\nHyrax on BN-254: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, commitment_size::<_, Hyrax<G1Affine254>>(num_vars));
    }

    println!("\nLigero on BN-254::Fr: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, commitment_size::<_, Ligero<Fr254>>(num_vars));
    }

    println!("\n---------------- Proof size ----------------");

    println!("\nHyrax on BLS12-381: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tser: {} B", proof_size::<_, Hyrax<G1Affine381>>(num_vars));
    }

    println!("\nLigero on BLS12-381::Fr: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, proof_size::<_, Ligero<Fr381>>(num_vars));
    }

    println!("\nHyrax on BN-254: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, proof_size::<_, Hyrax<G1Affine254>>(num_vars));
    }

    println!("\nLigero on BN-254::Fr: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!("\tnum_vars: {}, size: {} B", num_vars, proof_size::<_, Ligero<Fr254>>(num_vars));
    }

}
