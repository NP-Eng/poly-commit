use crate::{
    PCCommitment, PCCommitterKey, PCPreparedCommitment, PCRandomness, PCUniversalParams,
    PCVerifierKey,
};
use ark_crypto_primitives::crh::CRHScheme;
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, Path, TwoToOneParam};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;

/// The public parameters for any linear code PCS.
/// This is only a default setup with reasonable parameters.
/// To create your own public parameters, use:
/// # Example
/// ```rust
/// use ark_bls12_377::Fr;
/// use ark_crypto_primitives::{
///     crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
///     merkle_tree::{ByteDigestConverter, Config},
/// };
/// use ark_std::test_rng;
/// use ark_poly_commit::linear_codes::LinCodePCUniversalParams;
/// use core::marker::PhantomData;
///
/// type LeafH = Sha256;
/// type CompressH = Sha256;
/// struct MerkleTreeParams;
/// impl Config for MerkleTreeParams {
///     type Leaf = [u8];
///     type LeafDigest = <LeafH as CRHScheme>::Output;
///     type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
///     type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;
///     type LeafHash = LeafH;
///     type TwoToOneHash = CompressH;
/// }
/// type MTConfig = MerkleTreeParams;
/// let mut rng = &mut test_rng();
/// let leaf_hash_params = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
/// let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
///     .unwrap()
///     .clone();
/// let pp: LinCodePCUniversalParams<Fr, MTConfig> = LinCodePCUniversalParams::new(128, 2, true,
///     leaf_hash_params, two_to_one_params);
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCUniversalParams<F: PrimeField, C: Config, H: CRHScheme> {
    _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of the code rate.
    pub(crate) rho_inv: usize,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    #[derivative(Debug = "ignore")]
    pub(crate) col_hash_params: H::Parameters,
}

impl<F, C, H> LinCodePCUniversalParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    /// Create new UniversalParams
    pub fn new(
        sec_param: usize,
        rho_inv: usize,
        check_well_formedness: bool,
        leaf_hash_params: LeafParam<C>,
        two_to_one_params: TwoToOneParam<C>,
        col_hash_params: H::Parameters,
    ) -> Self {
        Self {
            _field: PhantomData,
            sec_param,
            rho_inv,
            check_well_formedness,
            leaf_hash_params,
            two_to_one_params,
            col_hash_params,
        }
    }
}

impl<F, C, H> PCUniversalParams for LinCodePCUniversalParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        if F::TWO_ADICITY < self.rho_inv as u32 {
            0
        } else if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }
}

/// Linear code commitment structure
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCCommitterKey<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    pub(crate) _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of code rate
    pub(crate) rho_inv: usize,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    /// Parameters for hashing the columns with a CRH Scheme
    #[derivative(Debug = "ignore")]
    pub(crate) col_hash_params: H::Parameters,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
}

impl<F, C, H> PCCommitterKey for LinCodePCCommitterKey<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        self.max_degree()
    }
}
/// The verifier key which holds some scheme parameters
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCVerifierKey<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    pub(crate) _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of code rate
    pub(crate) rho_inv: usize,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    /// Parameters for hashing the columns with a CRH Scheme
    #[derivative(Debug = "ignore")]
    pub(crate) col_hash_params: H::Parameters,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
}

impl<F, C, H> PCVerifierKey for LinCodePCVerifierKey<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        self.max_degree()
    }
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
/// Metadata about the PCS
pub struct Metadata {
    /// The number of rows of the matrix
    pub n_rows: usize,
    /// The number of columns of the matrix
    pub n_cols: usize,
    /// The number of columns of the matrix after encoding
    pub n_ext_cols: usize,
}

/// The commitment to a polynomial is a root of the merkle tree,
/// where each node is a hash of the column of the encoded coefficient matrix U.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCCommitment<C: Config> {
    /// number of rows resp. columns of the square matrix containing the coefficients of the polynomial
    pub metadata: Metadata,
    /// The actual commitment
    pub root: C::InnerDigest,
}

impl<C: Config> PCCommitment for LinCodePCCommitment<C> {
    fn empty() -> Self {
        LinCodePCCommitment::default()
    }

    fn has_degree_bound(&self) -> bool {
        false
    }
}

pub(crate) type LinCodePCPreparedCommitment<C> = LinCodePCCommitment<C>;

impl<Unprepared: PCCommitment, C: Config> PCPreparedCommitment<Unprepared>
    for LinCodePCPreparedCommitment<C>
{
    fn prepare(_cm: &Unprepared) -> Self {
        LinCodePCPreparedCommitment::default()
    }
}

pub(crate) type LinCodePCRandomness = ();

impl PCRandomness for LinCodePCRandomness {
    fn empty() -> Self {
        unimplemented!()
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Self {
        unimplemented!()
    }
}

/// Proof of an individual linear code well-formedness check or opening
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCProofSingle<F, C>
where
    F: PrimeField,
    C: Config,
{
    /// For each of the indices in q, `paths` contains the path from the root of the merkle tree to the leaf
    pub paths: Vec<Path<C>>,

    /// v, s.t. E(v) = w
    pub v: Vec<F>,

    /// Columns that are opened
    pub columns: Vec<Vec<F>>,
}

/// The Proof type for linear code PCS, which amounts to an array of individual proofs
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCProof<F, C>
where
    F: PrimeField,
    C: Config,
{
    /// Single opening
    pub opening: LinCodePCProofSingle<F, C>,
    pub(crate) well_formedness: Option<Vec<F>>,
}

// Multiple poly at one point
pub(crate) type LPCPArray<F, C> = Vec<LinCodePCProof<F, C>>;
