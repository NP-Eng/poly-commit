use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{
    merkle_tree::{Config, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::PrimeField;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_std::fmt::Debug;
use core::marker::PhantomData;
use digest::Digest;
use jf_primitives::pcs::transcript::IOPTranscript;
use std::borrow::Borrow;

use crate::ligero::utils::{inner_product, reed_solomon};
use crate::{Error, LabeledCommitment, LabeledPolynomial, PolynomialCommitment};

use ark_std::rand::RngCore;

mod utils;
use utils::Matrix;

mod data_structures;
use data_structures::*;

pub use data_structures::{Ligero, LigeroPCCommitterKey, LigeroPCProof, LigeroPCVerifierKey};

use utils::{calculate_t, compute_dimensions, get_indices_from_transcript, hash_column};

mod tests;

impl<F, C, D, S, P, const RHO_INV: usize, const SEC_PARAM: usize>
    Ligero<F, C, D, S, P, RHO_INV, SEC_PARAM>
where
    F: PrimeField,
    C: Config,
    Vec<u8>: Borrow<C::Leaf>,
    C::InnerDigest: Absorb,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
{
    /// Create a new instance of Ligero.
    /// If either or both parameters are None, their default values are used.
    pub fn new() -> Self {
        Self {
            _config: PhantomData,
            _field: PhantomData,
            // TODO potentially can get rid of digest and sponge
            _digest: PhantomData,
            _sponge: PhantomData,
            _poly: PhantomData,
        }
    }

    // /// The verifier can check the well-formedness of the commitment by taking random linear combinations.
    // fn check_well_formedness(
    //     root: &C::InnerDigest,
    //     well_formedness_proof: &LigeroPCProofSingle<F, C>,
    //     leaf_hash_params: &LeafParam<C>,
    //     two_to_one_params: &TwoToOneParam<C>,
    // ) -> Result<(), Error> {
    //     let t = calculate_t(RHO_INV, SEC_PARAM);

    //     // TODO replace unwraps by proper error handling
    //     let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"well_formedness_transcript");
    //     transcript
    //         .append_serializable_element(b"root", root)
    //         .unwrap();

    //     // 2. Get the linear combination coefficients from the transcript
    //     let mut r = Vec::new();
    //     for _ in 0..commitment.n_rows {
    //         r.push(transcript.get_and_append_challenge(b"r").unwrap());
    //     }
    //     // Upon sending `v` to the Verifier, add it to the sponge. Claim is that v = r.M
    //     transcript
    //         .append_serializable_element(b"v", &well_formedness_proof.v)
    //         .unwrap();

    //     Self::check_random_linear_combination(
    //         &r,
    //         &well_formedness_proof,
    //         &commitment.root,
    //         commitment.n_ext_cols,
    //         t,
    //         &mut transcript,
    //         leaf_hash_params,
    //         two_to_one_params,
    //     )
    // }

    fn check_random_linear_combination(
        coeffs: &[F],
        proof: &LigeroPCProofSingle<F, C>,
        root: &C::InnerDigest,
        n_ext_cols: usize,
        t: usize,
        transcript: &mut IOPTranscript<F>,
        leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> Result<(), Error> {
        // 1. Hash the received columns into leaf hashes
        let mut col_hashes: Vec<Vec<u8>> = Vec::new();
        for c in proof.columns.iter() {
            col_hashes.push(hash_column::<D, F>(c));
        }

        // 2. Compute t column indices to check the linear combination at
        let indices = get_indices_from_transcript::<F>(n_ext_cols, t, transcript);

        // 3. Verify the paths for each of the leaf hashes
        for (i, (leaf, q_i)) in col_hashes.into_iter().zip(indices.iter()).enumerate() {
            // TODO handle the error here
            let path = &proof.paths[i];
            assert!(
                path.leaf_index == *q_i,
                "Path is for a different index: i: {}, leaf index: {}!",
                q_i,
                path.leaf_index
            ); // TODO return an error

            path.verify(leaf_hash_params, two_to_one_params, root, leaf.clone())
                .unwrap();
        }

        // 4. Compute the encoding w = E(v)
        let w = reed_solomon(&proof.v, RHO_INV);

        // 5. Verify the random linear combinations
        for (transcript_index, matrix_index) in indices.into_iter().enumerate() {
            if inner_product(coeffs, &proof.columns[transcript_index]) != w[matrix_index] {
                // TODO return proper error
                return Err(Error::IncorrectInputLength(
                    "Incorrect linear combination".to_string(),
                ));
            }
        }

        Ok(())
    }
    fn compute_matrices(polynomial: &P) -> (Matrix<F>, Matrix<F>) {
        let mut coeffs = polynomial.coeffs().to_vec();

        // 1. Computing parameters and initial matrix
        let (n_rows, n_cols) = compute_dimensions::<F>(polynomial.degree() + 1); // for 6 coefficients, this is returning 4 x 2 with a row of 0s: fix

        // padding the coefficient vector with zeroes
        // TODO is this the most efficient/safest way to do it?
        coeffs.resize(n_rows * n_cols, F::zero());

        let mat = Matrix::new_from_flat(n_rows, n_cols, &coeffs);

        // 2. Apply Reed-Solomon encoding row-wise
        let ext_mat = Matrix::new_from_rows(
            mat.rows()
                .iter()
                .map(|r| reed_solomon(r, RHO_INV))
                .collect(),
        );

        (mat, ext_mat)
    }
    fn create_merkle_tree(
        ext_mat: &Matrix<F>,
        leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> MerkleTree<C> {
        let mut col_hashes: Vec<Vec<u8>> = Vec::new();
        let ext_mat_cols = ext_mat.cols();

        for col in ext_mat_cols.iter() {
            col_hashes.push(hash_column::<D, F>(col));
        }

        // pad the column hashes with zeroes
        let next_pow_of_two = col_hashes.len().next_power_of_two();
        col_hashes.resize(next_pow_of_two, vec![0; <D as Digest>::output_size()]);

        MerkleTree::<C>::new(leaf_hash_params, two_to_one_params, col_hashes).unwrap()
    }
    fn generate_proof(
        coeffs: &[F],
        mat: &Matrix<F>,
        ext_mat: &Matrix<F>,
        col_tree: &MerkleTree<C>,
        transcript: &mut IOPTranscript<F>,
    ) -> LigeroPCProofSingle<F, C> {
        let t = calculate_t(RHO_INV, SEC_PARAM, F::MODULUS_BIT_SIZE, ext_mat.m); // TODO this function will now probably need to take into account the number of rows/cols of the extended matrix

        // 1. Compute the linear combination using the random coefficients
        let v = mat.row_mul(coeffs);

        transcript.append_serializable_element(b"v", &v).unwrap();

        // 2. Generate t column indices to test the linear combination on
        let indices = get_indices_from_transcript(ext_mat.m, t, transcript);

        // 3. Compute Merkle tree paths for the columns
        let mut queried_columns = Vec::new();
        let mut paths = Vec::new();

        let ext_mat_cols = ext_mat.cols();

        for i in indices {
            queried_columns.push(ext_mat_cols[i].clone());
            paths.push(col_tree.generate_proof(i).unwrap());
        }

        LigeroPCProofSingle {
            paths,
            v,
            columns: queried_columns,
        }
    }
}

impl<F, P, S, C, D, const RHO_INV: usize, const SEC_PARAM: usize> PolynomialCommitment<F, P, S>
    for Ligero<F, C, D, S, P, RHO_INV, SEC_PARAM>
where
    F: PrimeField,
    P: DenseUVPolynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    Vec<u8>: Borrow<C::Leaf>,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
    C::InnerDigest: Absorb,
    D: Digest,
{
    type UniversalParams = LigeroPCUniversalParams;

    type CommitterKey = LigeroPCCommitterKey<C>;

    type VerifierKey = LigeroPCVerifierKey<C>;

    type PreparedVerifierKey = LigeroPCPreparedVerifierKey;

    type Commitment = LigeroPCCommitment<C>;

    type PreparedCommitment = LigeroPCPreparedCommitment;

    type Randomness = LigeroPCRandomness;

    type Proof = LPCPArray<F, C>;

    type BatchProof = Vec<Self::Proof>;

    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        assert!(
            RHO_INV >= 1,
            "RHO_INV is the inverse of the rate and must be at least 1"
        );
        // The domain will have size m * RHO_INV, but we already have the first m elements
        GeneralEvaluationDomain::<F>::compute_size_of_domain(max_degree * (RHO_INV - 1))
            .ok_or(Error::UnsupportedDegreeBound(max_degree))?;

        LigeroPCUniversalParams::default();
        Ok(())
    }

    fn trim(
        _pp: &Self::UniversalParams,
        _supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        todo!();
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        let mut commitments = Vec::new();

        for labeled_polynomial in polynomials.into_iter() {
            let polynomial = labeled_polynomial.polynomial();

            // 1. Compute matrices
            let (mat, ext_mat) = Self::compute_matrices(polynomial);

            // 2. Create the Merkle tree from the hashes of the columns
            let col_tree =
                Self::create_merkle_tree(&ext_mat, &ck.leaf_hash_params, &ck.two_to_one_params);

            // 3. Add root to transcript and generate random linear combination with it
            let root = col_tree.root();

            let mut transcript: IOPTranscript<F> =
                IOPTranscript::new(b"well_formedness_transcript");
            transcript
                .append_serializable_element(b"root", &root)
                .unwrap();

            let n_rows = mat.n;
            let n_cols = mat.m;
            let n_ext_cols = ext_mat.m;

            let mut r = Vec::new();
            for _ in 0..n_rows {
                r.push(transcript.get_and_append_challenge(b"r").unwrap());
            }

            // 4. Generate the proof by choosing random columns and proving their paths in the tree
            let well_formedness_proof = if ck.check_well_formedness {
                Some(Self::generate_proof(
                    &r,
                    &mat,
                    &ext_mat,
                    &col_tree,
                    &mut transcript,
                ))
            } else {
                None
            };

            let commitment = LigeroPCCommitment {
                n_rows,
                n_cols,
                n_ext_cols,
                root,
            };

            commitments.push(LabeledCommitment::new(
                labeled_polynomial.label().clone(),
                commitment,
                None, // TODO think about this (degree_bound)
            ));
        }

        Ok((commitments, Vec::new()))
        // TODO when should this return Err?
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        _challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        _rands: impl IntoIterator<Item = &'a Self::Randomness>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        let mut proof_array = LPCPArray::default();
        let labeled_commitments: Vec<&'a LabeledCommitment<Self::Commitment>> =
            commitments.into_iter().collect();
        let labeled_polynomials: Vec<&'a LabeledPolynomial<F, P>> =
            labeled_polynomials.into_iter().collect();

        assert_eq!(
            labeled_commitments.len(),
            labeled_polynomials.len(),
            // maybe return Err?
            "Mismatched lengths: {} commitments, {} polynomials",
            labeled_commitments.len(),
            labeled_polynomials.len()
        );

        for i in 0..labeled_polynomials.len() {
            let polynomial = labeled_polynomials[i].polynomial();
            let commitment = labeled_commitments[i].commitment();

            // TODO we receive a list of polynomials and a list of commitments
            // are we to understand that the first commitment is for the first polynomial, ...etc?

            // TODO we should maybe check that these two lists match, but that would imply recomputing merkle trees...
            // at least check labels?

            // 1. Compute matrices
            let (mat, ext_mat) = Self::compute_matrices(polynomial);

            // 2. Create the Merkle tree from the hashes of the columns
            let col_tree =
                Self::create_merkle_tree(&ext_mat, &ck.leaf_hash_params, &ck.two_to_one_params);

            // 3. Generate vector b and add v = b·M to the transcript
            let mut b = Vec::new();
            let point_pow = point.pow([commitment.n_cols as u64]); // TODO this and other conversions could potentially fail
            let mut acc_b = F::one();
            for _ in 0..commitment.n_rows {
                b.push(acc_b);
                acc_b *= point_pow;
            }

            let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"opening_transcript");

            transcript
                .append_serializable_element(b"point", point)
                .unwrap();

            proof_array.push(LigeroPCProof {
                opening: Self::generate_proof(&b, &mat, &ext_mat, &col_tree, &mut transcript),
                well_formedness: None,
            });
        }

        Ok(proof_array)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof_array: &Self::Proof,
        _challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let labeled_commitments: Vec<&'a LabeledCommitment<Self::Commitment>> =
            commitments.into_iter().collect();
        let values: Vec<F> = values.into_iter().collect();

        if labeled_commitments.len() != proof_array.len()
            || labeled_commitments.len() != values.len()
        {
            // maybe return Err?
            panic!("Mismatched lengths: {} proofs were provided for {} commitments with {} claimed values",
            labeled_commitments.len(), proof_array.len(), values.len());
        }

        for (i, labeled_commitment) in labeled_commitments.iter().enumerate() {
            let commitment = labeled_commitment.commitment();

            // TODO maybe check that the parameters have been calculated honestly (n_rows/cols/ext_cols);
            //      could they be used to cheat?
            let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"opening_transcript");

            // check if we've seen this commitment before. If not, we should verify it.
            if vk.check_well_formedness {
                if let None = &proof_array[i].well_formedness {
                    panic!("Handle the panic properly");
                }
                let well_formedness = &proof_array[i].well_formedness.as_ref().unwrap();
                transcript
                    .append_serializable_element(b"root", &commitment.root)
                    .unwrap();

                // 2. Get the linear combination coefficients from the transcript
                let mut r = Vec::new();
                for _ in 0..commitment.n_rows {
                    r.push(transcript.get_and_append_challenge(b"r").unwrap());
                }
                // Upon sending `v` to the Verifier, add it to the sponge. Claim is that v = r.M
                transcript
                    .append_serializable_element(b"v", &well_formedness.v)
                    .unwrap();

                // Self::check_random_linear_combination(
                //     &r,
                //     &well_formedness_proof,
                //     &commitment.root,
                //     commitment.n_ext_cols,
                //     t,
                //     &mut transcript,
                //     leaf_hash_params,
                //     two_to_one_params,
                // );
            }

            // 1. Compute a and b
            let mut a = Vec::new();
            let mut acc_a = F::one();
            for _ in 0..commitment.n_cols {
                a.push(acc_a);
                acc_a *= point;
            }

            // by now acc_a = point^n_cols
            let mut b = Vec::new();
            let mut acc_b = F::one();
            for _ in 0..commitment.n_rows {
                b.push(acc_b);
                acc_b *= acc_a;
            }
            let t = calculate_t(
                RHO_INV,
                SEC_PARAM,
                F::MODULUS_BIT_SIZE,
                commitment.n_ext_cols,
            ); // TODO include in ck/vk?

            // 2. Seed the transcript with the point and generate t random indices
            // TODO replace unwraps by proper error handling
            transcript
                .append_serializable_element(b"point", point)
                .unwrap();
            transcript
                .append_serializable_element(b"v", &proof_array[i].opening.v)
                .unwrap();

            // 3. Check the linear combination in the proof
            if {
                let coeffs: &[F] = &b;
                let root = &commitment.root;
                let n_ext_cols = commitment.n_ext_cols;
                let leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters = &vk.leaf_hash_params;
                let two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters = &vk.two_to_one_params;
                // 1. Hash the received columns into leaf hashes
                let mut col_hashes: Vec<Vec<u8>> = Vec::new();
                for c in proof_array[i].opening.columns.iter() {
                    col_hashes.push(hash_column::<D, F>(c));
                }

                // 2. Compute t column indices to check the linear combination at
                let indices = get_indices_from_transcript::<F>(n_ext_cols, t, &mut transcript);

                // 3. Verify the paths for each of the leaf hashes
                for (i, (leaf, q_i)) in col_hashes.into_iter().zip(indices.iter()).enumerate() {
                    // TODO handle the error here
                    let path = &proof_array[i].opening.paths[i];
                    assert!(
                        path.leaf_index == *q_i,
                        "Path is for a different index: i: {}, leaf index: {}!",
                        q_i,
                        path.leaf_index
                    ); // TODO return an error

                    path.verify(leaf_hash_params, two_to_one_params, root, leaf.clone())
                        .unwrap();
                }

                // 4. Compute the encoding w = E(v)
                let w = reed_solomon(&proof_array[i].opening.v, RHO_INV);


                // 5. Verify the random linear combinations
                for (transcript_index, matrix_index) in indices.into_iter().enumerate() {
                    if inner_product(coeffs, &proof_array[i].opening.columns[transcript_index]) != w[matrix_index] {
                        // TODO return proper error
                        return Err(Error::IncorrectInputLength(
                            "Incorrect linear combination".to_string(),
                        ));
                    }
                }

                Ok::<(), Self::Error>(())
            }
            .is_err()
            {
                // I think this can never be called since check_random_linear_combination will panick itself; must improve error handling
                println!("Function check failed verification of opening with index {i}");
                return Ok(false);
            }

            if inner_product(&proof_array[i].opening.v, &a) != values[i] {
                println!("Function check: claimed value in position {i} does not match the evaluation of the committed polynomial in the same position");
                return Ok(false);
            }
        }

        Ok(true)
    }
}

// TODO start considering degree bound
