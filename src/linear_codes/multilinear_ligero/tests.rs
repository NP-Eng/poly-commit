#[cfg(test)]
mod tests {

    use core::borrow::Borrow;
    use core::marker::PhantomData;

    use crate::linear_codes::LinearCodePCS;
    use crate::to_bytes;
    use crate::RngCore;
    use crate::{
        challenge::ChallengeGenerator,
        linear_codes::{utils::*, LinCodePCUniversalParams, MultilinearLigero},
        LabeledPolynomial, PolynomialCommitment,
    };
    use ark_bls12_377::Fq;
    use ark_bls12_377::Fr;
    use ark_bls12_381::Fr as Fr381;
    use ark_crypto_primitives::Error;
    use ark_crypto_primitives::{
        crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
        merkle_tree::{ByteDigestConverter, Config},
        sponge::poseidon::PoseidonSponge,
    };
    use ark_ff::{Field, PrimeField};
    use ark_poly::evaluations::multivariate::{MultilinearExtension, SparseMultilinearExtension};
    use ark_serialize::CanonicalSerialize;
    use ark_std::test_rng;
    use blake2::Blake2s256;
    use digest::Digest;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl pedersen::Window for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type LeafH = LeafIdentityHasher;
    type CompressH = Sha256;
    type ColHasher<F, D> = FieldToBytesColHasher<F, D>;

    struct FieldToBytesColHasher<F, D>
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

        fn setup<R: RngCore>(_rng: &mut R) -> Result<Self::Parameters, Error> {
            Ok(())
        }

        fn evaluate<T: Borrow<Self::Input>>(
            _parameters: &Self::Parameters,
            input: T,
        ) -> Result<Self::Output, Error> {
            let mut dig = D::new();
            dig.update(to_bytes!(input.borrow()).unwrap());
            Ok(dig.finalize().to_vec())
        }
    }

    struct LeafIdentityHasher;

    impl CRHScheme for LeafIdentityHasher {
        type Input = Vec<u8>;
        type Output = Vec<u8>;
        type Parameters = ();

        fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
            Ok(())
        }

        fn evaluate<T: Borrow<Self::Input>>(
            _: &Self::Parameters,
            input: T,
        ) -> Result<Self::Output, Error> {
            Ok(input.borrow().to_vec().into())
        }
    }

    struct MerkleTreeParams; //<F, D>(PhantomData<(F, D)>);

    impl Config for MerkleTreeParams {
        type Leaf = Vec<u8>;

        type LeafDigest = <LeafH as CRHScheme>::Output;
        type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafH;
        type TwoToOneHash = CompressH;
    }

    type MTConfig = MerkleTreeParams; //<F, Blake2s256>;
    type Sponge<F> = PoseidonSponge<F>;

    type LigeroPCS<F> = LinearCodePCS<
        MultilinearLigero<F, MTConfig, Blake2s256, Sponge<F>, SparseMultilinearExtension<F>>,
        F,
        SparseMultilinearExtension<F>,
        Sponge<F>,
        MTConfig,
        Blake2s256,
        ColHasher<F, Blake2s256>,
    >;

    fn rand_poly<Fr: PrimeField>(
        _: usize,
        num_vars: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> SparseMultilinearExtension<Fr> {
        match num_vars {
            Some(n) => SparseMultilinearExtension::rand(n, rng),
            None => unimplemented!(), // should not happen in ML case!
        }
    }

    fn constant_poly<Fr: PrimeField>(
        _: usize,
        num_vars: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> SparseMultilinearExtension<Fr> {
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        match num_vars {
            Some(n) => {
                let points = vec![(1, Fr::rand(rng))];
                SparseMultilinearExtension::from_evaluations(n, &points)
            }
            None => unimplemented!(), // should not happen in ML case!
        }
    }

    #[test]
    fn test_construction() {
        let mut rng = &mut test_rng();
        // just to make sure we have the right degree given the FFT domain for our field
        let leaf_hash_params = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();
        let col_hash_params = <ColHasher<Fr, Blake2s256> as CRHScheme>::setup(&mut rng).unwrap();
        let check_well_formedness = true;

        let pp: LinCodePCUniversalParams<Fr, MTConfig, ColHasher<Fr, Blake2s256>> =
            LinCodePCUniversalParams::new(
                128,
                4,
                check_well_formedness,
                leaf_hash_params,
                two_to_one_params,
                col_hash_params,
            );

        let (ck, vk) = LigeroPCS::<Fr>::trim(&pp, 0, 0, None).unwrap();

        let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
        let labeled_poly = LabeledPolynomial::new(
            "test".to_string(),
            rand_poly(1, Some(5), rand_chacha),
            Some(5),
            Some(5),
        );

        let mut test_sponge = test_sponge::<Fr>();
        let (c, rands) = LigeroPCS::<Fr>::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = rand_point(Some(5), rand_chacha);

        let value = labeled_poly.evaluate(&point);

        let mut challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
            ChallengeGenerator::new_univariate(&mut test_sponge);

        let proof = LigeroPCS::<Fr>::open(
            &ck,
            &[labeled_poly],
            &c,
            &point,
            &mut (challenge_generator.clone()),
            &rands,
            None,
        )
        .unwrap();
        assert!(LigeroPCS::<Fr>::check(
            &vk,
            &c,
            &point,
            [value],
            &proof,
            &mut challenge_generator,
            None
        )
        .unwrap());
    }

    #[test]
    fn test_calculate_t_with_good_parameters() {
        assert!(calculate_t::<Fq>(128, 4, 2_usize.pow(32)).unwrap() < 200);
        assert!(calculate_t::<Fq>(256, 4, 2_usize.pow(32)).unwrap() < 400);
    }

    #[test]
    fn test_calculate_t_with_bad_parameters() {
        calculate_t::<Fq>((Fq::MODULUS_BIT_SIZE - 60) as usize, 4, 2_usize.pow(60)).unwrap_err();
        calculate_t::<Fq>(400, 4, 2_usize.pow(32)).unwrap_err();
    }

    fn rand_point<F: Field>(num_vars: Option<usize>, rng: &mut ChaCha20Rng) -> Vec<F> {
        match num_vars {
            Some(n) => (0..n).map(|_| F::rand(rng)).collect(),
            None => unimplemented!(), // should not happen!
        }
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS<Fr>, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(10),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS<Fr>, _>(
            Some(10),
            constant_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(5),
            constant_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, LigeroPCS<Fr>, _>(
            Some(8),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(3),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, LigeroPCS<Fr>, _>(
            Some(10),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(5),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, LigeroPCS<Fr>, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(10),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, LigeroPCS<Fr>, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(8),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
