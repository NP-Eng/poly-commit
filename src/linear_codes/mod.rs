use crate::linear_codes::utils::*;
use crate::utils::ceil_div;
use crate::{
    Error, LabeledCommitment, LabeledPolynomial, PCCommitterKey, PCUniversalParams, PCVerifierKey,
    PolynomialCommitment,
};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::MerkleTree;
use ark_crypto_primitives::{merkle_tree::Config, sponge::CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial};
use ark_std::borrow::Borrow;
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::string::ToString;
use ark_std::vec::Vec;
#[cfg(not(feature = "std"))]
use num_traits::Float;

mod utils;

mod multilinear_ligero;
mod univariate_ligero;

pub use multilinear_ligero::MultilinearLigero;
pub use univariate_ligero::UnivariateLigero;

mod breakdown;
mod data_structures;
mod ligero;
use data_structures::*;

pub use data_structures::{BreakdownPCParams, LigeroPCParams, LinCodePCProof};

use utils::{calculate_t, get_indices_from_transcript, hash_column};

const FIELD_SIZE_ERROR: &str = "This field is not suitable for the proposed parameters";

/// This trait is another kir for this kiri interface
pub trait LinCodeInfo<C, H>
where
    C: Config,
    H: CRHScheme,
{
    /// Get security parameter
    fn sec_param(&self) -> usize;

    /// Get the inverse of code rate
    fn rho_inv(&self) -> (usize, usize);

    /// See whether there should be a well-formedness check
    fn check_well_formedness(&self) -> bool;

    /// Get LeafHash parameters
    fn leaf_hash_params(&self) -> &<<C as Config>::LeafHash as CRHScheme>::Parameters;

    /// Get TwoToOneHash paramters
    fn two_to_one_params(&self) -> &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters;

    /// Get column hashing parameters
    fn col_hash_params(&self) -> &H::Parameters;
}

/// A trait for linear encoding a messsage.
pub trait LinearEncode<F, C, P, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
    P: Polynomial<F>,
{
    /// For schemes like Breakdown and Ligero, PCCommiiterKey and
    /// PCVerifierKey and PCUniversalParams are all the same.
    type LinCodePCParams: PCUniversalParams + PCCommitterKey + PCVerifierKey + LinCodeInfo<C, H>;

    /// Does a default setup for the PCS.
    fn setup<R: RngCore>(
        rng: &mut R,
        leaf_hash_params: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
        col_hash_params: H::Parameters,
    ) -> Self::LinCodePCParams;

    /// Encode a message, which is interpreted as a vector of coefficients
    /// of a polynomial of degree m - 1.
    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Vec<F>;

    /// Get the representation of the polynomial
    fn poly_repr(polynomial: &P) -> Vec<F>;

    /// How we choose to split the query point into a Vec of Field elements.
    /// Needed for appending to transcript.
    fn point_to_vec(point: P::Point) -> Vec<F>;

    /// Compute the dimensions of an FFT-friendly (over F) matrix with at least n entries.
    /// The return pair (n, m) corresponds to the dimensions n x m.
    fn compute_dimensions(n: usize) -> (usize, usize) {
        assert_eq!(
            (n as f64) as usize,
            n,
            "n cannot be converted to f64: aborting"
        );

        let aux = (n as f64).sqrt().ceil() as usize;
        let n_cols = GeneralEvaluationDomain::<F>::new(aux)
            .expect("Field F does not admit FFT with m elements")
            .size();

        (ceil_div(n, n_cols), n_cols)
    }

    /// Compute the matrices for the polynomial
    fn compute_matrices(polynomial: &P, param: &Self::LinCodePCParams) -> (Matrix<F>, Matrix<F>) {
        let mut coeffs = Self::poly_repr(polynomial);

        // 1. Computing parameters and initial matrix
        let (n_rows, n_cols) = Self::compute_dimensions(coeffs.len()); // for 6 coefficients, this is returning 4 x 2 with a row of 0s: fix

        // padding the coefficient vector with zeroes
        // TODO is this the most efficient/safest way to do it?
        coeffs.resize(n_rows * n_cols, F::zero());

        let mat = Matrix::new_from_flat(n_rows, n_cols, &coeffs);

        // 2. Apply Reed-Solomon encoding row-wise
        let ext_mat =
            Matrix::new_from_rows(mat.rows().iter().map(|r| Self::encode(r, param)).collect());

        (mat, ext_mat)
    }

    /// Tensor the point
    fn tensor(point: &P::Point, left_len: usize, right_len: usize) -> (Vec<F>, Vec<F>);
}

/// Any linear-code-based commitment scheme.
pub struct LinearCodePCS<L, F, P, S, C, H>
where
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: Polynomial<F>,
    H: CRHScheme,
    L: LinearEncode<F, C, P, H>,
{
    _phantom: PhantomData<(L, F, P, S, C, H)>,
}

impl<L, F, P, S, C, H> PolynomialCommitment<F, P, S> for LinearCodePCS<L, F, P, S, C, H>
where
    L: LinearEncode<F, C, P, H>,
    F: PrimeField,
    P: Polynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    Vec<F>: Borrow<<H as CRHScheme>::Input>,
    H::Output: Into<C::Leaf>,
    C::Leaf: Sized + Clone + Default,
    H: CRHScheme,
{
    type UniversalParams = L::LinCodePCParams;

    type CommitterKey = L::LinCodePCParams;

    type VerifierKey = L::LinCodePCParams;

    type PreparedVerifierKey = LinCodePCPreparedVerifierKey;

    type Commitment = LinCodePCCommitment<C>;

    type PreparedCommitment = LinCodePCPreparedCommitment<C>;

    type Randomness = LinCodePCRandomness;

    type Proof = LPCPArray<F, C>;

    type BatchProof = Vec<Self::Proof>;

    type Error = Error;

    /// This is only a default setup with reasonable parameters.
    /// To create your own public parameters (from which vk/ck can be derived by `trim`),
    /// see the documentation for `LigeroPCUniversalParams`.
    fn setup<R: RngCore>(
        max_degree: usize,
        _num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let leaf_hash_params = <C::LeafHash as CRHScheme>::setup(rng).unwrap();
        let two_to_one_params = <C::TwoToOneHash as TwoToOneCRHScheme>::setup(rng)
            .unwrap()
            .clone();
        let col_hash_params = <H as CRHScheme>::setup(rng).unwrap();
        let pp = L::setup::<R>(rng, leaf_hash_params, two_to_one_params, col_hash_params);
        let real_max_degree = <Self::UniversalParams as PCUniversalParams>::max_degree(&pp);
        if max_degree > real_max_degree || real_max_degree == 0 {
            return Err(Error::InvalidParameters(FIELD_SIZE_ERROR.to_string()));
        }
        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        _supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        if <Self::UniversalParams as PCUniversalParams>::max_degree(pp) == 0 {
            return Err(Error::InvalidParameters(FIELD_SIZE_ERROR.to_string()));
        }
        Ok((pp.clone(), pp.clone()))
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

            // 1. Arrange the coefficients of the polynomial into a matrix,
            // and apply Reed-Solomon encoding to get `ext_mat`.
            let (mat, ext_mat) = L::compute_matrices(polynomial, ck);

            // 2. Create the Merkle tree from the hashes of each column.
            let col_tree = create_merkle_tree::<F, C, H>(
                &ext_mat,
                ck.leaf_hash_params(),
                ck.two_to_one_params(),
                ck.col_hash_params(),
            )?;

            // 3. Obtain the MT root and add it to the transcript.
            let root = col_tree.root();

            let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"transcript");

            transcript
                .append_serializable_element(b"root", &root)
                .map_err(|_| Error::TranscriptError)?;

            let n_rows = mat.n;
            let n_cols = mat.m;
            let n_ext_cols = ext_mat.m;

            // 4. The commitment is just the root, but since each commitment could be to a differently-sized polynomial, we also add some metadata.
            let commitment = LinCodePCCommitment {
                metadata: Metadata {
                    n_rows,
                    n_cols,
                    n_ext_cols,
                },
                root,
            };

            commitments.push(LabeledCommitment::new(
                labeled_polynomial.label().clone(),
                commitment,
                None,
            ));
        }
        let com_len = &commitments.len();
        Ok((commitments, vec![(); *com_len]))
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

        if labeled_commitments.len() != labeled_polynomials.len() {
            return Err(Error::IncorrectInputLength(format!(
                "Mismatched lengths: {} commitments, {} polynomials",
                labeled_commitments.len(),
                labeled_polynomials.len()
            )));
        }

        for i in 0..labeled_polynomials.len() {
            let polynomial = labeled_polynomials[i].polynomial();
            let commitment = labeled_commitments[i].commitment();
            let n_rows = commitment.metadata.n_rows;
            let n_cols = commitment.metadata.n_cols;
            let root = &commitment.root;

            // 1. Arrange the coefficients of the polynomial into a matrix,
            // and apply Reed-Solomon encoding to get `ext_mat`.
            let (mat, ext_mat) = L::compute_matrices(polynomial, ck);

            // 2. Create the Merkle tree from the hashes of each column.
            let col_tree = create_merkle_tree::<F, C, H>(
                &ext_mat,
                ck.leaf_hash_params(),
                ck.two_to_one_params(),
                ck.col_hash_params(),
            )?;

            // 3. Generate vector `b = [1, z^m, z^(2m), ..., z^((m-1)m)]`
            // This could potentially fail when n_cols > 1<<64, but `ck` won't allow commiting to such polynomials.
            // let point_pow = point.pow([n_cols as u64]);
            let (_, b) = L::tensor(point, n_cols, n_rows);

            let mut transcript = IOPTranscript::new(b"transcript");
            transcript
                .append_serializable_element(b"root", root)
                .map_err(|_| Error::TranscriptError)?;

            // If we are checking well-formedness, we need to compute the well-formedness proof (which is just r.M) and append it to the transcript.
            let well_formedness = if ck.check_well_formedness() {
                let mut r = Vec::new();
                for _ in 0..n_rows {
                    r.push(
                        transcript
                            .get_and_append_challenge(b"r")
                            .map_err(|_| Error::TranscriptError)?,
                    );
                }
                let v = mat.row_mul(&r);

                transcript
                    .append_serializable_element(b"v", &v)
                    .map_err(|_| Error::TranscriptError)?;
                Some(v)
            } else {
                None
            };

            let point_vec = L::point_to_vec(point.clone());
            for element in point_vec.iter() {
                transcript
                    .append_serializable_element(b"point", element)
                    .map_err(|_| Error::TranscriptError)?;
            }

            proof_array.push(LinCodePCProof {
                // compute the opening proof and append b.M to the transcript
                opening: generate_proof(
                    ck.sec_param(),
                    ck.rho_inv(),
                    &b,
                    &mat,
                    &ext_mat,
                    &col_tree,
                    &mut transcript,
                )?,
                well_formedness,
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
            return Err(Error::IncorrectInputLength(
                format!(
                    "Mismatched lengths: {} proofs were provided for {} commitments with {} claimed values",labeled_commitments.len(), proof_array.len(), values.len()
                )
            ));
        }
        let leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters =
            vk.leaf_hash_params();
        let two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters =
            vk.two_to_one_params();

        for (i, labeled_commitment) in labeled_commitments.iter().enumerate() {
            let commitment = labeled_commitment.commitment();
            let n_rows = commitment.metadata.n_rows;
            let n_cols = commitment.metadata.n_cols;
            let n_ext_cols = commitment.metadata.n_ext_cols;
            let root = &commitment.root;
            let t = calculate_t::<F>(vk.sec_param(), vk.rho_inv(), n_ext_cols)?;

            let mut transcript = IOPTranscript::new(b"transcript");
            transcript
                .append_serializable_element(b"root", &commitment.root)
                .map_err(|_| Error::TranscriptError)?;

            let out = if vk.check_well_formedness() {
                if proof_array[i].well_formedness.is_none() {
                    return Err(Error::InvalidCommitment);
                }
                let tmp = &proof_array[i].well_formedness.as_ref();
                let well_formedness = tmp.unwrap();
                let mut r = Vec::with_capacity(n_rows);
                for _ in 0..n_rows {
                    r.push(
                        transcript
                            .get_and_append_challenge(b"r")
                            .map_err(|_| Error::TranscriptError)?,
                    );
                }
                // Upon sending `v` to the Verifier, add it to the sponge. Claim is that v = r.M
                transcript
                    .append_serializable_element(b"v", well_formedness)
                    .map_err(|_| Error::TranscriptError)?;

                (Some(well_formedness), Some(r))
            } else {
                (None, None)
            };

            // 1. Seed the transcript with the point and the recieved vector
            // TODO Consider removing the evaluation point from the transcript.
            let point_vec = L::point_to_vec(point.clone());
            for element in point_vec.iter() {
                transcript
                    .append_serializable_element(b"point", element)
                    .map_err(|_| Error::TranscriptError)?;
            }
            transcript
                .append_serializable_element(b"v", &proof_array[i].opening.v)
                .map_err(|_| Error::TranscriptError)?;

            // 2. Ask random oracle for the `t` indices where the checks happen
            let indices = get_indices_from_transcript::<F>(n_ext_cols, t, &mut transcript)?;

            // 3. Hash the received columns into leaf hashes
            let col_hashes: Vec<C::Leaf> = proof_array[i]
                .opening
                .columns
                .iter()
                .map(|c| hash_column::<F, C, H>(c.clone(), &vk.col_hash_params()).unwrap())
                .collect();

            // 4. Verify the paths for each of the leaf hashes - this is only run once,
            // even if we have a well-formedness check (i.e., we save sending and checking the columns).
            // See "Concrete optimizations to the commitment scheme", p.12 of [Brakedown](https://eprint.iacr.org/2021/1043.pdf)
            for (j, (leaf, q_j)) in col_hashes.iter().zip(indices.iter()).enumerate() {
                let path = &proof_array[i].opening.paths[j];
                if path.leaf_index != *q_j {
                    return Err(Error::InvalidCommitment);
                }

                path.verify(leaf_hash_params, two_to_one_params, root, leaf.clone())
                    .map_err(|_| Error::InvalidCommitment)?;
            }

            // helper closure: checks if a.b = c
            let check_inner_product = |a, b, c| -> Result<(), Error> {
                if inner_product(a, b) != c {
                    return Err(Error::InvalidCommitment);
                }

                Ok(())
            };

            // 5. Compute the encoding w = E(v)
            let w = L::encode(&proof_array[i].opening.v, vk);

            // 6. Compute a = [1, z, z^2, ..., z^(n_cols_1)]
            // where z denotes the query `point`.
            let (a, b) = L::tensor(point, n_cols, n_rows);

            // Compute b = [1, z^n_cols, z^(2*n_cols), ..., z^((n_rows-1)*n_cols)]
            let coeffs: &[F] = &b;

            // 7. Probabilistic checks that whatever the prover sent,
            // matches with what the verifier computed for himself.
            // Note: we sacrifice some code repetition in order not to repeat execution.
            if let (Some(well_formedness), Some(r)) = out {
                let w_well_formedness = L::encode(well_formedness, vk);
                for (transcript_index, matrix_index) in indices.iter().enumerate() {
                    check_inner_product(
                        &r,
                        &proof_array[i].opening.columns[transcript_index],
                        w_well_formedness[*matrix_index],
                    )?;
                    check_inner_product(
                        coeffs,
                        &proof_array[i].opening.columns[transcript_index],
                        w[*matrix_index],
                    )?;
                }
            } else {
                for (transcript_index, matrix_index) in indices.iter().enumerate() {
                    check_inner_product(
                        coeffs,
                        &proof_array[i].opening.columns[transcript_index],
                        w[*matrix_index],
                    )?;
                }
            }

            if inner_product(&proof_array[i].opening.v, &a) != values[i] {
                eprintln!("Function check: claimed value in position {i} does not match the evaluation of the committed polynomial in the same position");
                return Ok(false);
            }
        }

        Ok(true)
    }
}

// TODO maybe this can go to utils
fn create_merkle_tree<F, C, H>(
    ext_mat: &Matrix<F>,
    leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
    two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    col_hash_params: &H::Parameters,
) -> Result<MerkleTree<C>, Error>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
    Vec<F>: Borrow<<H as CRHScheme>::Input>,
    H::Output: Into<C::Leaf>,
    C::Leaf: Default + Clone,
{
    let mut col_hashes: Vec<C::Leaf> = Vec::new();
    let ext_mat_cols = ext_mat.cols();

    for col in ext_mat_cols.into_iter() {
        let col_digest = hash_column::<F, C, H>(col, &col_hash_params)?;
        col_hashes.push(col_digest);
    }

    // pad the column hashes with zeroes
    let next_pow_of_two = col_hashes.len().next_power_of_two();
    col_hashes.resize(next_pow_of_two, <C::Leaf>::default());

    MerkleTree::<C>::new(leaf_hash_params, two_to_one_params, col_hashes)
        .map_err(|_| Error::HashingError)
}

fn generate_proof<F, C>(
    sec_param: usize,
    rho_inv: (usize, usize),
    b: &[F],
    mat: &Matrix<F>,
    ext_mat: &Matrix<F>,
    col_tree: &MerkleTree<C>,
    transcript: &mut IOPTranscript<F>,
) -> Result<LinCodePCProofSingle<F, C>, Error>
where
    F: PrimeField,
    C: Config,
{
    let t = calculate_t::<F>(sec_param, rho_inv, ext_mat.n)?;

    // 1. left-multiply the matrix by `b`, where for a requested query point `z`,
    // `b = [1, z^m, z^(2m), ..., z^((m-1)m)]`
    let v = mat.row_mul(b);

    transcript
        .append_serializable_element(b"v", &v)
        .map_err(|_| Error::TranscriptError)?;

    // 2. Generate t column indices to test the linear combination on
    let indices = get_indices_from_transcript(ext_mat.m, t, transcript)?;

    // 3. Compute Merkle tree paths for the requested columns
    let mut queried_columns = Vec::with_capacity(t);
    let mut paths = Vec::with_capacity(t);

    let ext_mat_cols = ext_mat.cols();

    for i in indices {
        queried_columns.push(ext_mat_cols[i].clone());
        paths.push(
            col_tree
                .generate_proof(i)
                .map_err(|_| Error::TranscriptError)?,
        );
    }

    Ok(LinCodePCProofSingle {
        paths,
        v,
        columns: queried_columns,
    })
}
