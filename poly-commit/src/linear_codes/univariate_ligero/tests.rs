#[cfg(test)]
mod tests {
    use crate::{
        linear_codes::LigeroPCParams,
        test_types::{
            test_sponge, FieldToBytesColHasher, LeafIdentityHasher, TestMerkleTreeParams,
            TestUVLigero,
        },
        LabeledPolynomial, PolynomialCommitment,
    };

    use ark_bls12_377::Fr;
    use ark_bls12_381::Fr as Fr381;
    use ark_crypto_primitives::crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme};
    use ark_ff::{Field, PrimeField, UniformRand};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
    use ark_std::test_rng;
    use blake2::Blake2s256;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

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
        let mut rng = &mut test_rng();
        // just to make sure we have the right degree given the FFT domain for our field
        let leaf_hash_param = <LeafIdentityHasher as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_hash_param = <Sha256 as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();
        let col_hash_params =
            <FieldToBytesColHasher<Fr, Blake2s256> as CRHScheme>::setup(&mut rng).unwrap();
        let check_well_formedness = true;

        let pp: LigeroPCParams<Fr, TestMerkleTreeParams, FieldToBytesColHasher<Fr, Blake2s256>> =
            LigeroPCParams::new(
                128,
                4,
                check_well_formedness,
                leaf_hash_param,
                two_to_one_hash_param,
                col_hash_params,
            );

        let (ck, vk) = TestUVLigero::<Fr>::trim(&pp, 0, 0, None).unwrap();

        let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
        let labeled_poly = LabeledPolynomial::new(
            "test".to_string(),
            rand_poly(degree, None, rand_chacha),
            None,
            None,
        );

        let mut test_sponge = test_sponge::<Fr>();
        let (c, rands) = TestUVLigero::<Fr>::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = Fr::rand(rand_chacha);

        let value = labeled_poly.evaluate(&point);

        let proof = TestUVLigero::<Fr>::open(
            &ck,
            &[labeled_poly],
            &c,
            &point,
            &mut (test_sponge.clone()),
            &rands,
            None,
        )
        .unwrap();
        assert!(TestUVLigero::<Fr>::check(
            &vk,
            &c,
            &point,
            [value],
            &proof,
            &mut test_sponge,
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
        single_poly_test::<_, _, TestUVLigero<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, TestUVLigero<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, TestUVLigero<Fr>, _>(
            None,
            constant_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, TestUVLigero<Fr381>, _>(
            None,
            constant_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, TestUVLigero<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, TestUVLigero<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, _, TestUVLigero<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        linear_poly_degree_bound_test::<_, _, TestUVLigero<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, _, TestUVLigero<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, _, TestUVLigero<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, _, TestUVLigero<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, _, TestUVLigero<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, _, TestUVLigero<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, _, TestUVLigero<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, TestUVLigero<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, TestUVLigero<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, TestUVLigero<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, TestUVLigero<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, TestUVLigero<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, TestUVLigero<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, _, TestUVLigero<Fr>, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_degree_bound_test::<_, _, TestUVLigero<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, TestUVLigero<Fr>, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test::<Fr>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, TestUVLigero<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test::<Fr381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    #[should_panic]
    fn bad_degree_bound_test() {
        use crate::tests::*;
        use ark_bls12_381::Fq as Fq381;
        bad_degree_bound_test::<_, _, TestUVLigero<Fq381>, _>(
            rand_poly::<Fq381>,
            rand_point::<Fq381>,
            poseidon_sponge_for_test::<Fq381>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
    }
}
