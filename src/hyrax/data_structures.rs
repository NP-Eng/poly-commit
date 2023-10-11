

use ark_ec::AffineRepr;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use ark_std::rand::RngCore;

use crate::{PCUniversalParams, PCCommitterKey, PCPreparedVerifierKey, PCPreparedCommitment, PCRandomness};

/// TODO should the length be contained in any of these structures?

/// `UniversalParams` amount to a Pederson commitment key
/// TODO do we want this derivative business?
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct HyraxUniversalParams<G: AffineRepr> {
    /// A list of generators of the group.
    pub com_key: Vec<G>,

    /// A generator of the group.
    pub h: G,
}

impl<G: AffineRepr> PCUniversalParams for HyraxUniversalParams<G> {
    fn max_degree(&self) -> usize {
        // Only MLEs are supported
        1
    }
}

/// The committer key is used to commit to scalars, and by extension, to commit
/// to polynomials and open those commitments
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Default(bound = ""),
    // TODO remove or re-introduce
    // Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
// Cannot simply alias HyraxUniversalParams because of the derivation
pub struct HyraxCommitterKey<G: AffineRepr> {
    /// A list of generators of the group.
    pub com_key: Vec<G>,

    /// A generator of the group.
    pub h: G,
}

/// The verifier key, which coincides with the committer key
pub type HyraxVerifierKey<G> = HyraxCommitterKey<G>;

impl<G: AffineRepr> PCCommitterKey for HyraxCommitterKey<G> {
    fn max_degree(&self) -> usize {
        // Only MLEs are supported
        1
    }
    fn supported_degree(&self) -> usize {
        // Only MLEs are supported
        1
    }
}

/// Nothing to do to prepare this prover-verifier key.
pub type HyraxPreparedVerifierKey<G> = HyraxVerifierKey<G>;

impl<G: AffineRepr> PCPreparedVerifierKey<HyraxVerifierKey<G>> for HyraxPreparedVerifierKey<G> {
    /// Simply clone the prover-verifier key
    fn prepare(vk: &HyraxVerifierKey<G>) -> Self {
        vk.clone()
    }
}

pub struct HyraxCommitment<G: AffineRepr> {
    /// A list of multi-commits to each row of the matrix containing the
    /// polynomial.
    pub row_coms: Vec<G>,
}

pub type HyraxPreparedCommitment<E> = HyraxCommitment<E>;

impl<G: AffineRepr> PCPreparedCommitment<HyraxCommitment<G>> for HyraxPreparedCommitment<G> {
    /// Simply clone the prover-verifier key
    fn prepare(vk: &HyraxCommitment<G>) -> Self {
        vk.clone()
    }
}

pub(crate) type HyraxRandomness = ();

/// This object is not used in the Hyrax PCS (instead, the Pedersen
/// commitment module generates the necessary randomness).
impl PCRandomness for HyraxRandomness {
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

pub struct HyraxProof<G: AffineRepr> {
    com_eval: G,
    com_d: G,
    com_b: G,
}
