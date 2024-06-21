use crate::{
    linear_codes::{LinearCodePCS, MultilinearBrakedown, MultilinearLigero, UnivariateLigero},
    to_bytes,
};
use ark_crypto_primitives::{
    crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, SparseMultilinearExtension};
use ark_serialize::CanonicalSerialize;
use ark_std::{borrow::Borrow, marker::PhantomData, rand::RngCore};
use blake2::Blake2s256;
use digest::Digest;

type LeafH = LeafIdentityHasher;
type CompressH = Sha256;
type ColHasher<F, D> = FieldToBytesColHasher<F, D>;

pub struct LeafIdentityHasher;

impl CRHScheme for LeafIdentityHasher {
    type Input = Vec<u8>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        Ok(input.borrow().to_vec().into())
    }
}

/// Needed for benches and tests.
pub struct FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    _phantom: PhantomData<(F, D)>,
}

impl<F, D> CRHScheme for FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    type Input = Vec<F>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_rng: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _parameters: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        let mut dig = D::new();
        dig.update(to_bytes!(input.borrow()).unwrap());
        Ok(dig.finalize().to_vec())
    }
}

pub struct TestMerkleTreeParams;

impl Config for TestMerkleTreeParams {
    type Leaf = Vec<u8>;

    type LeafDigest = <LeafH as CRHScheme>::Output;
    type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

pub type TestUVLigero<F> = LinearCodePCS<
    UnivariateLigero<F, TestMerkleTreeParams, DensePolynomial<F>, ColHasher<F, Blake2s256>>,
    F,
    DensePolynomial<F>,
    TestMerkleTreeParams,
    ColHasher<F, Blake2s256>,
>;

pub type TestMLLigero<F> = LinearCodePCS<
    MultilinearLigero<
        F,
        TestMerkleTreeParams,
        SparseMultilinearExtension<F>,
        ColHasher<F, Blake2s256>,
    >,
    F,
    SparseMultilinearExtension<F>,
    TestMerkleTreeParams,
    ColHasher<F, Blake2s256>,
>;

pub type TestMLBrakedown<F> = LinearCodePCS<
    MultilinearBrakedown<
        F,
        TestMerkleTreeParams,
        SparseMultilinearExtension<F>,
        ColHasher<F, Blake2s256>,
    >,
    F,
    SparseMultilinearExtension<F>,
    TestMerkleTreeParams,
    ColHasher<F, Blake2s256>,
>;
