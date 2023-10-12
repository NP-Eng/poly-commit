#[cfg(test)]
mod tests {

    use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
    use ark_ec::AffineRepr;
    use ark_ed_on_bls12_381::EdwardsAffine;
    use ark_ff::PrimeField;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    use crate::challenge::ChallengeGenerator;
    use crate::hyrax::HyraxPC;

    use crate::linear_codes::utils::test_sponge;
    use crate::{LabeledPolynomial, PolynomialCommitment};

    #[inline]
    fn rand_dense_multilinear<F: PrimeField>(
        num_vars: usize,
        rng: &mut ChaCha20Rng,
    ) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::rand(num_vars, rng)
    }

    #[inline]
    fn rand_point<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
        (0..num_vars).map(|_| F::rand(rng)).collect()
    }

    #[test]
    fn test_hyrax_construction() {
        // Desired number of variables (must be even!)
        let n = 8;

        type Hyrax = HyraxPC<EdwardsAffine>;
        type Fr = <EdwardsAffine as AffineRepr>::ScalarField;

        let chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

        let pp = Hyrax::setup(1, Some(8), chacha).unwrap();

        let (ck, vk) = Hyrax::trim(&pp, 1, 1, None).unwrap();

        let l_poly = LabeledPolynomial::new(
            "test_poly".to_string(),
            rand_dense_multilinear::<Fr>(n, chacha),
            None,
            None,
        );

        let (c, rands) = Hyrax::commit(&ck, &[l_poly.clone()], Some(chacha)).unwrap();

        let point: Vec<Fr> = rand_point(n, chacha);
        let value = l_poly.evaluate(&point);

        // Dummy argument
        let mut test_sponge = test_sponge::<Fr>();
        let mut challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
            ChallengeGenerator::new_univariate(&mut test_sponge);

        let proof = Hyrax::open(
            &ck,
            &[l_poly],
            &c,
            &point,
            &mut (challenge_generator.clone()),
            &rands,
            Some(chacha),
        )
        .unwrap();

        assert!(Hyrax::check(
            &vk,
            &c,
            &point,
            [value],
            &proof,
            &mut challenge_generator,
            Some(chacha),
        )
        .unwrap());
    }
}
