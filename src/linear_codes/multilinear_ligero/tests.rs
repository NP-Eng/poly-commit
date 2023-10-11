#[cfg(test)]
mod tests {

    use core::borrow::Borrow;
    use core::marker::PhantomData;

    use crate::linear_codes::LinearCodePCS;
    use crate::tests::poseidon_parameters_for_test;
    use crate::RngCore;
    use crate::{
        challenge::ChallengeGenerator,
        linear_codes::{utils::*, LinCodePCUniversalParams, MultilinearLigero},
        LabeledPolynomial, PolynomialCommitment,
    };
    use ark_bn254::Fr;
    use ark_crypto_primitives::merkle_tree::IdentityDigestConverter;
    use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
    use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
    use ark_crypto_primitives::Error;
    use ark_crypto_primitives::{
        crh::{CRHScheme, TwoToOneCRHScheme},
        merkle_tree::Config,
        sponge::poseidon::PoseidonSponge,
    };
    use ark_ff::{Field, PrimeField};
    use ark_poly::evaluations::multivariate::{MultilinearExtension, SparseMultilinearExtension};
    use ark_std::test_rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    // We introduce the wrapper only for the purpose of `setup` function not panicing with unimplemented
    struct PoseidonWrapper<F>(PhantomData<F>);

    impl<F: PrimeField + Absorb> CRHScheme for PoseidonWrapper<F> {
        type Input = [F];
        type Output = F;
        type Parameters = PoseidonConfig<F>;

        fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
            Ok(poseidon_parameters_for_test())
        }

        fn evaluate<T: Borrow<Self::Input>>(
            parameters: &Self::Parameters,
            input: T,
        ) -> Result<Self::Output, Error> {
            let input = input.borrow();

            let mut sponge = PoseidonSponge::new(parameters);
            sponge.absorb(&input);
            let res = sponge.squeeze_field_elements::<F>(1);
            Ok(res[0])
        }
    }

    impl<F: PrimeField + Absorb> TwoToOneCRHScheme for PoseidonWrapper<F> {
        type Input = F;
        type Output = F;
        type Parameters = PoseidonConfig<F>;

        fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
            Ok(poseidon_parameters_for_test())
        }

        fn evaluate<T: Borrow<Self::Input>>(
            parameters: &Self::Parameters,
            left_input: T,
            right_input: T,
        ) -> Result<Self::Output, Error> {
            let left_input = left_input.borrow();
            let right_input = right_input.borrow();

            let mut sponge = PoseidonSponge::new(parameters);
            sponge.absorb(left_input);
            sponge.absorb(right_input);
            let res = sponge.squeeze_field_elements::<F>(1);
            Ok(res[0])
        }

        fn compress<T: Borrow<Self::Output>>(
            parameters: &Self::Parameters,
            left_input: T,
            right_input: T,
        ) -> Result<Self::Output, Error> {
            let left_input = left_input.borrow();
            let right_input = right_input.borrow();

            let mut sponge = PoseidonSponge::new(parameters);
            sponge.absorb(left_input);
            sponge.absorb(right_input);
            let res = sponge.squeeze_field_elements::<F>(1);
            Ok(res[0])
        }
    }

    type ColHasher<F> = PoseidonWrapper<F>;
    type CompressH<F> = PoseidonWrapper<F>;

    struct LeafIdentityHasher<F>(PhantomData<F>);

    impl<F: PrimeField> CRHScheme for LeafIdentityHasher<F> {
        type Input = F;
        type Output = F;
        type Parameters = ();

        fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, Error> {
            Ok(())
        }

        fn evaluate<T: Borrow<Self::Input>>(
            _: &Self::Parameters,
            input: T,
        ) -> Result<Self::Output, Error> {
            Ok(*input.borrow())
        }
    }

    struct MerkleTreeParams<F>(PhantomData<F>);

    impl<F: PrimeField + Absorb> Config for MerkleTreeParams<F> {
        type Leaf = F;

        type LeafDigest = <LeafIdentityHasher<F> as CRHScheme>::Output;
        type LeafInnerDigestConverter = IdentityDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH<F> as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafIdentityHasher<F>;
        type TwoToOneHash = CompressH<F>;
    }

    type MTConfig<F> = MerkleTreeParams<F>; 
    type Sponge<F> = PoseidonSponge<F>;

    type LigeroPCS<F> = LinearCodePCS<
        MultilinearLigero<F, MTConfig<F>, Sponge<F>, SparseMultilinearExtension<F>>,
        F,
        SparseMultilinearExtension<F>,
        Sponge<F>,
        MTConfig<F>,
        ColHasher<F>,
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
        // just to make sure we have the right degree given the FFT domain for our field
        let leaf_hash_params = ();
        let col_hash_params = poseidon_parameters_for_test();
        let two_to_one_params = poseidon_parameters_for_test();
        let check_well_formedness = true;

        let pp: LinCodePCUniversalParams<Fr, MTConfig<Fr>, ColHasher<Fr>> =
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
    }
}
