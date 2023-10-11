#[cfg(test)]
mod tests {

    use core::borrow::Borrow;
    use core::marker::PhantomData;

    use crate::ark_std::UniformRand;
    use crate::linear_codes::LinearCodePCS;
    use crate::PolynomialCommitment;
    use crate::{
        challenge::ChallengeGenerator,
        linear_codes::{utils::*, LinCodePCUniversalParams, UnivariateLigero},
        LabeledPolynomial,
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
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
    use ark_std::test_rng;
    use rand_chacha::rand_core::RngCore;
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
        UnivariateLigero<F, MTConfig<F>, Sponge<F>, DensePolynomial<F>>,
        F,
        DensePolynomial<F>,
        Sponge<F>,
        MTConfig<F>,
        ColHasher<F>,
    >;

    fn rand_poly<Fr: PrimeField>(
        degree: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<Fr> {
        DensePolynomial::rand(degree, rng)
    }

    fn constant_poly<Fr: PrimeField>(
        _: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<Fr> {
        DensePolynomial::from_coefficients_slice(&[Fr::rand(rng)])
    }

    #[test]
    fn test_construction() {
        let degree = 4;
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
            rand_poly(degree, None, rand_chacha),
            None,
            None,
        );

        let mut test_sponge = test_sponge::<Fr>();
        let (c, rands) = LigeroPCS::<Fr>::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = Fr::rand(rand_chacha);

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

    fn rand_point<F: Field>(_: Option<usize>, rng: &mut ChaCha20Rng) -> F {
        F::rand(rng)
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS<Fr>, _>(
            None,
            constant_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, LigeroPCS<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, _, LigeroPCS<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, _, LigeroPCS<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, _, LigeroPCS<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, _, LigeroPCS<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, LigeroPCS<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, LigeroPCS<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, LigeroPCS<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, _, LigeroPCS<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, LigeroPCS<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bn254");
    }
}
