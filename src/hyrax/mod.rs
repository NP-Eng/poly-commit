
mod data_structures;
mod utils;
pub use data_structures::*;

use core::marker::PhantomData;
use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, CryptographicSponge};
use ark_ec::{CurveGroup, AffineRepr, VariableBaseMSM};
use ark_poly::{DenseMultilinearExtension, Polynomial};
use ark_std::rand::RngCore;
use blake2::Blake2s256;

use crate::{PolynomialCommitment, Error, LabeledPolynomial, LabeledCommitment, hyrax::utils::flat_to_matrix_column_major, challenge::ChallengeGenerator};

/// Hyrax polynomial committment scheme:
/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups. This is a
/// Fiat-Shamired version of the PCS described in the Hyrax paper
/// [[WTsTW17]][hyrax], with the difference that, unlike in the
/// cited reference, the evaluation of the polynomial at the point
/// of interest is indeed revealed to the verifier at the end.
///
/// [hyrax]: https://eprint.iacr.org/2017/1132.pdf
pub struct HyraxPC<
    // The curve used for Pedersen commitments (only EC groups are
    // supported as of now).
    G: AffineRepr,
    // TODO make generic or fix one type and remove this
    // S: CryptographicSponge,
> {
    _curve: PhantomData<G>,
}


// TODO ********************************************************

impl<G: AffineRepr> HyraxPC<G> {
    fn pedersen_commit(
        com_key: &HyraxCommitterKey<G>,
        scalars: &[G::ScalarField],
    ) -> (G, G::ScalarField) {
        // This block is taken from the function `cm_commit` in the IPA
        // module from this crate
        let scalars_bigint = ark_std::cfg_iter!(scalars)
            .map(|s| s.into_bigint())
            .collect::<Vec<_>>();
        let mut com = <G::Group as VariableBaseMSM>::msm_bigint(com_key, &scalars_bigint);

        let r = G::ScalarField::rand();
        com += &com_key.h.mul(r);

        (com, r)
    }
} 

type MLE<G: AffineRepr> = DenseMultilinearExtension<G::ScalarField>;

// TODO ********************************************************

impl<G:AffineRepr> PolynomialCommitment<
    G::ScalarField,
    DenseMultilinearExtension<G::ScalarField>,
    // Dummy sponge - required by the trait, not used in this implementation
    PoseidonSponge<G::ScalarField>,
> for HyraxPC<G> {
    type UniversalParams = HyraxUniversalParams<G>;
    type CommitterKey = HyraxCommitterKey<G>;
    type VerifierKey = HyraxVerifierKey<G>;
    type PreparedVerifierKey = HyraxPreparedVerifierKey<G>;
    type Commitment = HyraxCommitment<G>;
    type PreparedCommitment = HyraxPreparedCommitment<G>;
    type Randomness = HyraxRandomness<G>;
    type Proof = HyraxProof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Outputs mock universal parameters for the Hyrax polynomial commitment
    /// scheme. It does *not* return random keys across calls and should never
    /// be used in settings where security is required - it is only useful for
    /// testing. Furthermore, the point at infinity could possibly be part of
    /// the output, which sould not happen in an actual key.
    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {

        G::BaseField;
        
        let n = num_vars.expect("Hyrax requires num_vars to be specified");
        
        assert_eq!(max_degree, 1, "Hyrax only supports multilinear polynomials");
        assert_eq!(n % 2, 0, "Only polynomials with an even number of variables \
                    are supported in this implementation");

        // Number of rows (or, equivalently, colums) of a square matrix
        // containing the coefficients of an n-variate ML polynomial
        let dim = 1 << n / 2;

        // The following block of code is largely taking from the IPA module
        // in this crate.
        let points: Vec<_> = ark_std::cfg_into_iter!(0..dim + 1)
        .map(|i| {
            let i = i as u64;
            let mut hash =
                Blake2s256::digest(["Hyrax protocol", &i.to_le_bytes()].concat().as_slice());
            let mut p = G::from_random_bytes(&hash);
            let mut j = 0u64;
            while p.is_none() {
                // PROTOCOL NAME, i, j
                let mut bytes = "Hyrax protocol".to_vec();
                bytes.extend(i.to_le_bytes());
                bytes.extend(j.to_le_bytes());
                hash = Blake2s256::digest(bytes.as_slice());
                p = G::from_random_bytes(&hash);
                j += 1;
            }
            let point = p.unwrap();
            point.mul_by_cofactor_to_group()
        })
        .collect();

        G::Group::normalize_batch(&points);

        let h = points.pop().unwrap();

        Ok(HyraxUniversalParams { com_key: points, h, num_vars: n })
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        
        assert!(supported_degree == 1 && supported_hiding_bound = 1, 
            "Hyrax only supports multilinear polynomials: \
            the passed degrees should be 1");

        assert!(enforced_degree_bounds.is_none(),
            "Hyrax only supports multilinear polynomials: \
            enforced_degree_bounds should be `None`"); 
        
        let HyraxUniversalParams { com_key, h, num_vars } = pp;

        let ck = HyraxCommitterKey { com_key, h, num_vars };

        let vk: HyraxVerifierKey = ck.clone();

        Ok((ck, vk))
    }

    /// Outputs a list of commitments to the passed polynomials
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, MLE<G>>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        MLE<G>: 'a,
    {
        let mut coms = Vec::new();
        let mut rands = Vec::new();

        for l_poly in polynomials {

            let mut com_rands = Vec::new();

            let label = l_poly.label();
            let poly = l_poly.polynomial();

            assert_eq!(l_poly.degree_bound().unwrap_or(1), 1,
                "Hyrax only supports ML polynomials: the degree bound should \
                be Some(1) or None"
            );

            assert!(l_poly.hiding_bound().is_none(),
                "Hiding bounds are not part of the Hyrax PCS");

            let n = poly.num_vars();
            let dim = 1 << n / 2;

            assert!(
                n <= ck.num_vars,
                "Attempted to commit to a polynomial with {n} variables, but
                this key only supports up to {} variables",
                ck.num_vars
            );

            let m = flat_to_matrix_column_major(poly.to_evaluations(), dim, dim);

            let row_coms = m.map(|row| {
                let (c, r) = Self::pedersen_commit(ck, &row);
                com_rands.push(r);
                c
            });

            let com = HyraxCommitment { row_coms };
            let l_comm = LabeledCommitment::new(label.to_string(), com, 1);

            coms.push(l_comm);
            rands.push(com_rands);
        }

        Ok((coms, rands))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, MLE<G>>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a <DenseMultilinearExtension<G::ScalarField> as Polynomial<G::ScalarField>>::Point,
        // Not used and not generic on the cryptographic sponge S
        _opening_challenges: &mut ChallengeGenerator<G::ScalarField, PoseidonSponge<G::ScalarField>>,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::Randomness: 'a,
        MLE<G>: 'a,
    {
        let mut combined_polynomial = P::zero();
        let mut combined_rand = G::ScalarField::zero();
        let mut combined_commitment_proj = G::Group::zero();

        let mut has_hiding = false;

        let polys_iter = labeled_polynomials.into_iter();
        let rands_iter = rands.into_iter();
        let comms_iter = commitments.into_iter();

        let combine_time = start_timer!(|| "Combining polynomials, randomness, and commitments.");

        let mut cur_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);

        for (labeled_polynomial, (labeled_commitment, randomness)) in
            polys_iter.zip(comms_iter.zip(rands_iter))
        {
            let label = labeled_polynomial.label();
            assert_eq!(labeled_polynomial.label(), labeled_commitment.label());
            Self::check_degrees_and_bounds(ck.supported_degree(), labeled_polynomial)?;

            let polynomial = labeled_polynomial.polynomial();
            let degree_bound = labeled_polynomial.degree_bound();
            let hiding_bound = labeled_polynomial.hiding_bound();
            let commitment = labeled_commitment.commitment();

            combined_polynomial += (cur_challenge, polynomial);
            combined_commitment_proj += &commitment.comm.mul(cur_challenge);

            if hiding_bound.is_some() {
                has_hiding = true;
                combined_rand += &(cur_challenge * &randomness.rand);
            }

            cur_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);

            let has_degree_bound = degree_bound.is_some();

            assert_eq!(
                has_degree_bound,
                commitment.shifted_comm.is_some(),
                "shifted_comm mismatch for {}",
                label
            );

            assert_eq!(
                degree_bound,
                labeled_commitment.degree_bound(),
                "labeled_comm degree bound mismatch for {}",
                label
            );
            if let Some(degree_bound) = degree_bound {
                let shifted_polynomial = Self::shift_polynomial(ck, polynomial, degree_bound);
                combined_polynomial += (cur_challenge, &shifted_polynomial);
                combined_commitment_proj += &commitment.shifted_comm.unwrap().mul(cur_challenge);

                if hiding_bound.is_some() {
                    let shifted_rand = randomness.shifted_rand;
                    assert!(
                        shifted_rand.is_some(),
                        "shifted_rand.is_none() for {}",
                        label
                    );
                    combined_rand += &(cur_challenge * &shifted_rand.unwrap());
                }
            }

            cur_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);
        }

        end_timer!(combine_time);

        let combined_v = combined_polynomial.evaluate(point);

        // Pad the coefficients to the appropriate vector size
        let d = ck.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        let mut combined_commitment;
        let mut hiding_commitment = None;

        if has_hiding {
            let mut rng = rng.expect("hiding commitments require randomness");
            let hiding_time = start_timer!(|| "Applying hiding.");
            let mut hiding_polynomial = P::rand(d, &mut rng);
            hiding_polynomial -= &P::from_coefficients_slice(&[hiding_polynomial.evaluate(point)]);
            let hiding_rand = G::ScalarField::rand(&mut rng);
            let hiding_commitment_proj = Self::cm_commit(
                ck.comm_key.as_slice(),
                hiding_polynomial.coeffs(),
                Some(ck.s),
                Some(hiding_rand),
            );

            let mut batch =
                G::Group::normalize_batch(&[combined_commitment_proj, hiding_commitment_proj]);
            hiding_commitment = Some(batch.pop().unwrap());
            combined_commitment = batch.pop().unwrap();

            let mut byte_vec = Vec::new();
            combined_commitment
                .serialize_uncompressed(&mut byte_vec)
                .unwrap();
            point.serialize_uncompressed(&mut byte_vec).unwrap();
            combined_v.serialize_uncompressed(&mut byte_vec).unwrap();
            hiding_commitment
                .unwrap()
                .serialize_uncompressed(&mut byte_vec)
                .unwrap();
            let bytes = byte_vec.as_slice();
            let hiding_challenge = Self::compute_random_oracle_challenge(bytes);
            combined_polynomial += (hiding_challenge, &hiding_polynomial);
            combined_rand += &(hiding_challenge * &hiding_rand);
            combined_commitment_proj +=
                &(hiding_commitment.unwrap().mul(hiding_challenge) - &ck.s.mul(combined_rand));

            end_timer!(hiding_time);
        }

        let combined_rand = if has_hiding {
            Some(combined_rand)
        } else {
            None
        };

        let proof_time =
            start_timer!(|| format!("Generating proof for degree {} combined polynomial", d + 1));

        combined_commitment = combined_commitment_proj.into_affine();

        // ith challenge
        let mut byte_vec = Vec::new();
        combined_commitment
            .serialize_uncompressed(&mut byte_vec)
            .unwrap();
        point.serialize_uncompressed(&mut byte_vec).unwrap();
        combined_v.serialize_uncompressed(&mut byte_vec).unwrap();
        let bytes = byte_vec.as_slice();
        let mut round_challenge = Self::compute_random_oracle_challenge(bytes);

        let h_prime = ck.h.mul(round_challenge).into_affine();

        // Pads the coefficients with zeroes to get the number of coeff to be d+1
        let mut coeffs = combined_polynomial.coeffs().to_vec();
        if coeffs.len() < d + 1 {
            for _ in coeffs.len()..(d + 1) {
                coeffs.push(G::ScalarField::zero());
            }
        }
        let mut coeffs = coeffs.as_mut_slice();

        // Powers of z
        let mut z: Vec<G::ScalarField> = Vec::with_capacity(d + 1);
        let mut cur_z: G::ScalarField = G::ScalarField::one();
        for _ in 0..(d + 1) {
            z.push(cur_z);
            cur_z *= point;
        }
        let mut z = z.as_mut_slice();

        // This will be used for transforming the key in each step
        let mut key_proj: Vec<G::Group> = ck.comm_key.iter().map(|x| (*x).into()).collect();
        let mut key_proj = key_proj.as_mut_slice();

        let mut temp;

        // Key for MSM
        // We initialize this to capacity 0 initially because we want to use the key slice first
        let mut comm_key = &ck.comm_key;

        let mut l_vec = Vec::with_capacity(log_d);
        let mut r_vec = Vec::with_capacity(log_d);

        let mut n = d + 1;
        while n > 1 {
            let (coeffs_l, coeffs_r) = coeffs.split_at_mut(n / 2);
            let (z_l, z_r) = z.split_at_mut(n / 2);
            let (key_l, key_r) = comm_key.split_at(n / 2);
            let (key_proj_l, _) = key_proj.split_at_mut(n / 2);

            let l = Self::cm_commit(key_l, coeffs_r, None, None)
                + &h_prime.mul(Self::inner_product(coeffs_r, z_l));

            let r = Self::cm_commit(key_r, coeffs_l, None, None)
                + &h_prime.mul(Self::inner_product(coeffs_l, z_r));

            let lr = G::Group::normalize_batch(&[l, r]);
            l_vec.push(lr[0]);
            r_vec.push(lr[1]);

            let mut byte_vec = Vec::new();
            round_challenge
                .serialize_uncompressed(&mut byte_vec)
                .unwrap();
            lr[0].serialize_uncompressed(&mut byte_vec).unwrap();
            lr[1].serialize_uncompressed(&mut byte_vec).unwrap();
            let bytes = byte_vec.as_slice();
            round_challenge = Self::compute_random_oracle_challenge(bytes);
            let round_challenge_inv = round_challenge.inverse().unwrap();

            ark_std::cfg_iter_mut!(coeffs_l)
                .zip(coeffs_r)
                .for_each(|(c_l, c_r)| *c_l += &(round_challenge_inv * &*c_r));

            ark_std::cfg_iter_mut!(z_l)
                .zip(z_r)
                .for_each(|(z_l, z_r)| *z_l += &(round_challenge * &*z_r));

            ark_std::cfg_iter_mut!(key_proj_l)
                .zip(key_r)
                .for_each(|(k_l, k_r)| *k_l += &(k_r.mul(round_challenge)));

            coeffs = coeffs_l;
            z = z_l;

            key_proj = key_proj_l;
            temp = G::Group::normalize_batch(key_proj);
            comm_key = &temp;

            n /= 2;
        }

        end_timer!(proof_time);

        Ok(Proof {
            l_vec,
            r_vec,
            final_comm_key: comm_key[0],
            c: coeffs[0],
            hiding_comm: hiding_commitment,
            rand: combined_rand,
        })
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        opening_challenges: &mut ChallengeGenerator<G::ScalarField, S>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        let d = vk.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        if proof.l_vec.len() != proof.r_vec.len() || proof.l_vec.len() != log_d {
            return Err(Error::IncorrectInputLength(
                format!(
                    "Expected proof vectors to be {:}. Instead, l_vec size is {:} and r_vec size is {:}",
                    log_d,
                    proof.l_vec.len(),
                    proof.r_vec.len()
                )
            ));
        }

        let check_poly =
            Self::succinct_check(vk, commitments, *point, values, proof, opening_challenges);

        if check_poly.is_none() {
            return Ok(false);
        }

        let check_poly_coeffs = check_poly.unwrap().compute_coeffs();
        let final_key = Self::cm_commit(
            vk.comm_key.as_slice(),
            check_poly_coeffs.as_slice(),
            None,
            None,
        );
        if !(final_key - &proof.final_comm_key.into()).is_zero() {
            return Ok(false);
        }

        end_timer!(check_time);
        Ok(true)
    }

    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        values: &Evaluations<G::ScalarField, P::Point>,
        proof: &Self::BatchProof,
        opening_challenges: &mut ChallengeGenerator<G::ScalarField, S>,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, (point_label, point)) in query_set.iter() {
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }

        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut randomizer = G::ScalarField::one();

        let mut combined_check_poly = P::zero();
        let mut combined_final_key = G::Group::zero();

        for ((_point_label, (point, labels)), p) in query_to_labels_map.into_iter().zip(proof) {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", labels.len()));
            let mut comms: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut vals = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = values
                    .get(&(label.clone(), *point))
                    .ok_or(Error::MissingEvaluation {
                        label: label.to_string(),
                    })?;

                comms.push(commitment);
                vals.push(*v_i);
            }

            let check_poly = Self::succinct_check(
                vk,
                comms.into_iter(),
                *point,
                vals.into_iter(),
                p,
                opening_challenges,
            );

            if check_poly.is_none() {
                return Ok(false);
            }

            let check_poly = P::from_coefficients_vec(check_poly.unwrap().compute_coeffs());
            combined_check_poly += (randomizer, &check_poly);
            combined_final_key += &p.final_comm_key.mul(randomizer);

            randomizer = u128::rand(rng).into();
            end_timer!(lc_time);
        }

        let proof_time = start_timer!(|| "Checking batched proof");
        let final_key = Self::cm_commit(
            vk.comm_key.as_slice(),
            combined_check_poly.coeffs(),
            None,
            None,
        );
        if !(final_key - &combined_final_key).is_zero() {
            return Ok(false);
        }

        end_timer!(proof_time);

        Ok(true)
    }

    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<G::ScalarField>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &mut ChallengeGenerator<G::ScalarField, S>,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<G::ScalarField, Self::BatchProof>, Self::Error>
    where
        Self::Randomness: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        let label_poly_map = polynomials
            .into_iter()
            .zip(rands)
            .zip(commitments)
            .map(|((p, r), c)| (p.label(), (p, r, c)))
            .collect::<BTreeMap<_, _>>();

        let mut lc_polynomials = Vec::new();
        let mut lc_randomness = Vec::new();
        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();

        for lc in linear_combinations {
            let lc_label = lc.label().clone();
            let mut poly = P::zero();
            let mut degree_bound = None;
            let mut hiding_bound = None;

            let mut combined_comm = G::Group::zero();
            let mut combined_shifted_comm: Option<G::Group> = None;

            let mut combined_rand = G::ScalarField::zero();
            let mut combined_shifted_rand: Option<G::ScalarField> = None;

            let num_polys = lc.len();
            for (coeff, label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                let label: &String = label.try_into().expect("cannot be one!");
                let &(cur_poly, cur_rand, cur_comm) =
                    label_poly_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                if num_polys == 1 && cur_poly.degree_bound().is_some() {
                    assert!(
                        coeff.is_one(),
                        "Coefficient must be one for degree-bounded equations"
                    );
                    degree_bound = cur_poly.degree_bound();
                } else if cur_poly.degree_bound().is_some() {
                    eprintln!("Degree bound when number of equations is non-zero");
                    return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                }

                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());

                combined_rand += &(cur_rand.rand * coeff);
                combined_shifted_rand = Self::combine_shifted_rand(
                    combined_shifted_rand,
                    cur_rand.shifted_rand,
                    *coeff,
                );

                let commitment = cur_comm.commitment();
                combined_comm += &commitment.comm.mul(*coeff);
                combined_shifted_comm = Self::combine_shifted_comm(
                    combined_shifted_comm,
                    commitment.shifted_comm,
                    *coeff,
                );
            }

            let lc_poly =
                LabeledPolynomial::new(lc_label.clone(), poly, degree_bound, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_randomness.push(Randomness {
                rand: combined_rand,
                shifted_rand: combined_shifted_rand,
            });

            lc_commitments.push(combined_comm);
            if let Some(combined_shifted_comm) = combined_shifted_comm {
                lc_commitments.push(combined_shifted_comm);
            }

            lc_info.push((lc_label, degree_bound));
        }

        let lc_commitments = Self::construct_labeled_commitments(&lc_info, &lc_commitments);

        let proof = Self::batch_open(
            ck,
            lc_polynomials.iter(),
            lc_commitments.iter(),
            &query_set,
            opening_challenges,
            lc_randomness.iter(),
            rng,
        )?;
        Ok(BatchLCProof { proof, evals: None })
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<G::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        eqn_query_set: &QuerySet<P::Point>,
        eqn_evaluations: &Evaluations<P::Point, G::ScalarField>,
        proof: &BatchLCProof<G::ScalarField, Self::BatchProof>,
        opening_challenges: &mut ChallengeGenerator<G::ScalarField, S>,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let BatchLCProof { proof, .. } = proof;
        let label_comm_map = commitments
            .into_iter()
            .map(|c| (c.label(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();
        let mut evaluations = eqn_evaluations.clone();
        for lc in linear_combinations {
            let lc_label = lc.label().clone();
            let num_polys = lc.len();

            let mut degree_bound = None;
            let mut combined_comm = G::Group::zero();
            let mut combined_shifted_comm: Option<G::Group> = None;

            for (coeff, label) in lc.iter() {
                if label.is_one() {
                    for (&(ref label, _), ref mut eval) in evaluations.iter_mut() {
                        if label == &lc_label {
                            **eval -= coeff;
                        }
                    }
                } else {
                    let label: &String = label.try_into().unwrap();
                    let &cur_comm = label_comm_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                    if num_polys == 1 && cur_comm.degree_bound().is_some() {
                        assert!(
                            coeff.is_one(),
                            "Coefficient must be one for degree-bounded equations"
                        );
                        degree_bound = cur_comm.degree_bound();
                    } else if cur_comm.degree_bound().is_some() {
                        return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                    }

                    let commitment = cur_comm.commitment();
                    combined_comm += &commitment.comm.mul(*coeff);
                    combined_shifted_comm = Self::combine_shifted_comm(
                        combined_shifted_comm,
                        commitment.shifted_comm,
                        *coeff,
                    );
                }
            }

            lc_commitments.push(combined_comm);

            if let Some(combined_shifted_comm) = combined_shifted_comm {
                lc_commitments.push(combined_shifted_comm);
            }

            lc_info.push((lc_label, degree_bound));
        }

        let lc_commitments = Self::construct_labeled_commitments(&lc_info, &lc_commitments);

        Self::batch_check(
            vk,
            &lc_commitments,
            &eqn_query_set,
            &evaluations,
            proof,
            opening_challenges,
            rng,
        )
    }
}
