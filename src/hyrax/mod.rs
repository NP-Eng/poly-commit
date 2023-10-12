
mod data_structures;
mod utils;
pub use data_structures::*;

use core::marker::PhantomData;
use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, CryptographicSponge};
use ark_ec::{CurveGroup, AffineRepr, VariableBaseMSM};
use ark_poly::{DenseMultilinearExtension, Polynomial};
use ark_std::rand::RngCore;

use crate::linear_codes::utils::{Matrix, inner_product, vector_sum, scalar_by_vector};

use blake2::Blake2s256;

use crate::{PolynomialCommitment, Error, LabeledPolynomial, LabeledCommitment, hyrax::utils::{flat_to_matrix_column_major, usize_to_bits, naive_chi}, challenge::ChallengeGenerator};

/// Hyrax polynomial committment scheme:
/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups. This is a
/// Fiat-Shamired version of the PCS described in the Hyrax paper
/// [[WTsTW17]][hyrax].
///
/// [hyrax]: https://eprint.iacr.org/2017/1132.pdf
/// 
/// * Modification note *
/// 
/// In the PCS contained in the cited article, the verifier never learns the
/// actual evaluation of the polynomial at the requested point, but is instead
/// convinced that a previously received Pedersen commitment is indeed a
/// commitment to said evaluation - this is what the SNARK proposed therein
/// necessitates. However, the Arkworks framework requies the verifier to
/// actually learn that value, which is why we have added the opening of
/// the commitment at the end of the protocol. This might not result in an
/// optimal non-hiding PCS, but we feel it is the most faithful adaptation of
/// original PCS that can be implemented with the current restrictions.
pub struct HyraxPC<
    // The curve used for Pedersen commitments (only EC groups are
    // supported as of now).
    G: AffineRepr,
    // TODO make generic or fix one type and remove this
    // S: CryptographicSponge,
> {
    _curve: PhantomData<G>,
}

// TODO so far everything is done with asserts instead of the Error
// types defined by the library. Is this okay?

// TODO use ark_std::cfg_iter! instead of iter() as it is now?

// TODO add "trusting" version of utils linear-algebra functions which do not check dimensions?

// TODO ********************************************************

impl<G: AffineRepr> HyraxPC<G> {
    fn pedersen_commit(
        com_key: &HyraxCommitterKey<G>,
        scalars: &[G::ScalarField],
        r: Option<G::ScalarField>,
    ) -> (G, G::ScalarField) {
        // TODO does this trim the key suitably?
        // This block is taken from the function `cm_commit` in the IPA
        // module from this crate
        let scalars_bigint = ark_std::cfg_iter!(scalars)
            .map(|s| s.into_bigint())
            .collect::<Vec<_>>();
        let mut com = <G::Group as VariableBaseMSM>::msm_bigint(com_key, &scalars_bigint);

        let r = r.unwrap_or(G::ScalarField::rand());
        com += &com_key.h.mul(r);

        (com, r)
    }
} 

type MLE<G: AffineRepr> = DenseMultilinearExtension<G::ScalarField>;

// TODO ********************************************************

impl<G: AffineRepr> PolynomialCommitment<
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
    type Proof = Vec<HyraxProof<G>>;
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
                let (c, r) = Self::pedersen_commit(ck, &row, None);
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
        // TODO is it safe to open several polynomials at once?
        // TODO is there a more efficient way to open several polynomials at once?
        //      can one e.g. share zs, ds..?

        let n = point.len();

        assert_eq!(n % 2, 0, "Only points with an even number of variables \
            are supported in this implementation");

        let dim = 1 << n / 2;

        let point_lower = point[n / 2..];
        let point_upper = point[..n / 2];

        // TODO this way to compute the bits is very inefficient
        let l = (0..dim).map(|idx| naive_chi(&usize_to_bits(idx, n / 2), &point_lower));
        let r = (0..dim).map(|idx| naive_chi(&usize_to_bits(idx, n / 2), &point_upper));

        let proofs = Vec::new();

        for (l_poly, (l_com, randomness)) in
            labeled_polynomials.into_iter()
                .zip(commitments.into_iter()
                    .zip(rands.into_iter()))
        {

            // TODO check if the poly was actually necessary
            let label = l_poly.label();
            assert_eq!(label, l_com.label(), "Mismatching labels: {label} and {}", l_com.label());

            let poly = l_poly.polynomial();
            let com = l_com.commitment();

            // TODO chech num of vars matches n
            assert_eq!(poly.num_vars(), n, "The committed polynomial has {} variables, but \
                the point has {n} variables", poly.num_vars());

            let t_aux = flat_to_matrix_column_major(poly.to_evaluations(), dim, dim);
            let t = Matrix::new_from_rows(t_aux);

            let lt = t.row_mul(&l);

            let r_lt = l.iter().zip(r.iter()).map(|(l, r)| l * r).sum();

            // TODO change to multi-exponentiation OR directly compute as the commitment to LT?
            let t_prime = com.com_key.iter().zip(l.iter()).map(|(c, s)| c.mul(l)).sum();

            let eval = inner_product(&lt, r);

            // Singleton commit
            let (com_eval, r_eval) = Self::pedersen_commit(ck, &[eval], None);

            // ******** Dot product argument ********
            // Appendix A.2 in the reference article

            let rng = rng.expect("Opening polynomials requires randomness");

            let mut d: Vec<G::ScalarField> = (0..dim).map(|_| G::BaseField::rand(rng)).collect();

            let b = inner_product(&r, &d);

            // Multi-commit
            let (com_d, r_d) = Self::pedersen_commit(ck, &d, None);

            // Singleton commit
            let (com_b, r_b) = Self::pedersen_commit(ck, &[b], None);

            // Receive the random challenge c from the verifier, i.e. squeeze
            // it from the transcript.
            // TODO
            let c = G::BaseField::rand(rng);

            let z = vector_sum(&d, &scalar_by_vector(c, &lt));
            let z_d = c * r_lt + r_d;
            let z_b = c * r_eval + r_b;

            // ******** Opening ********
            // This is *not* part of the Hyrax PCS as described in the reference
            // article. Cf. the "Modification note" at the beginning of this file.
            // From the prover's perspective, opening amounts to adding r_eval to
            // the proof.

            {proofs.push(HyraxProof {com_eval, com_d, com_b, z, z_d, z_b, r_eval});            
        }

        Ok(proofs)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        _opening_challenges: &mut ChallengeGenerator<G::ScalarField, S>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {

        let n = point.len();

        assert_eq!(n % 2, 0, "Only points with an even number of variables \
            are supported in this implementation");

        let dim = 1 << n / 2;

        let point_lower = point[n / 2..];
        let point_upper = point[..n / 2];

        // TODO this way to compute the bits is very inefficient
        let l = (0..dim).map(|idx| naive_chi(&usize_to_bits(idx, n / 2), &point_lower));
        let r = (0..dim).map(|idx| naive_chi(&usize_to_bits(idx, n / 2), &point_upper));

        for (com, (claim, h_proof)) in
            commitments.into_iter()
                .zip(values.into_iter()
                    .zip(proof.into_iter()))
        {

            let row_coms = com.row_coms;

            // extract each field from h_proof
            let HyraxProof {com_eval, com_d, com_b, z, z_d, z_b, r_eval} = h_proof;

            // TODO chech num of vars matches n
            assert_eq!(row_coms.len(), 1 << n / 2, "The commitment should have 2^(n/2) = has {} entries, but \
                it has {} instead", 1 << n / 2, row_coms.len());

            // TODO change to multi-exponentiation OR directly compute as the commitment to LT?
            let t_prime = row_coms.iter().zip(l.iter()).map(|(c, s)| c.mul(l)).sum();

            // Receive the random challenge c from the verifier, i.e. squeeze
            // it from the transcript.
            // TODO
            let c = G::BaseField::rand(rng);

            // First check
            let (com_z_zd, _) = Self::pedersen_commit(vk, &z, z_d);
            if com_z_zd != c * t_prime + com_d {
                return Ok(false)
            }

            // Second check
            let (com_dp, _) = Self::pedersen_commit(vk, &inner_product(&r, &z), z_b);
            if com_dp != c * com_eval + com_b {
                return Ok(false)
            }

            // Third check: opening
            if com_eval != Self::pedersen_commit(vk, &[claim], r_eval) {
                return Ok(false)
            }
        }

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
