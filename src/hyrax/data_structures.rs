

use ark_ec::AffineRepr;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use ark_std::{rand::RngCore, UniformRand};



use crate::{PCUniversalParams, PCCommitterKey, PCPreparedVerifierKey, PCPreparedCommitment, PCRandomness, PCCommitment, PCVerifierKey};

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
    /// Maximum number of variables a polynomial can be committed to with this
    /// key
    pub num_vars: usize,
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
    /// Maximum number of variables a polynomial can be committed to with this
    /// key
    pub num_vars: usize,
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

impl<G: AffineRepr> PCVerifierKey for HyraxVerifierKey<G> {
    fn max_degree(&self) -> usize {
        1
    }

    fn supported_degree(&self) -> usize {
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

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct HyraxCommitment<G: AffineRepr> {
    /// A list of multi-commits to each row of the matrix containing the
    /// polynomial.
    pub row_coms: Vec<G>,
}

impl<G: AffineRepr> PCCommitment for HyraxCommitment<G> {
    #[inline]
    fn empty() -> Self {
        HyraxCommitment {
            row_coms: Vec::new(),
        }
    }

    // The degree bound is always 1, since only multilinear polynomials are
    // supported
    fn has_degree_bound(&self) -> bool {
        true
    }
}

pub type HyraxPreparedCommitment<E> = HyraxCommitment<E>;

impl<G: AffineRepr> PCPreparedCommitment<HyraxCommitment<G>> for HyraxPreparedCommitment<G> {
    /// Simply clone the prover-verifier key
    fn prepare(vk: &HyraxCommitment<G>) -> Self {
        vk.clone()
    }
}

pub(crate) type HyraxRandomness<G: AffineRepr> = Vec<G::ScalarField>;

/// A vector of scalars, each of which multiplies the distinguished group
/// element in the Pederson commitment key for a different commitment
impl<G: AffineRepr> PCRandomness for HyraxRandomness<G> {
    fn empty() -> Self {
        unimplemented!()
    }

    fn rand<R: RngCore>(
        num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        rng: &mut R,
    ) -> Self {
        (0..num_queries).map(|_| G::ScalarField::rand(rng)).collect()
    }
}

pub struct HyraxProof<G: AffineRepr> {
    pub com_eval: G,
    pub com_d: G,
    pub com_b: G,
    pub z: Vec<G::ScalarField>,
    pub z_d: G::ScalarField,
    pub z_b: G::ScalarField,
    // The seed r_eval is not part of a Hyrax PCS proof as described in the
    // reference article. Cf. the "Modification note" at the beginning of
    // mod.rs
    pub r_eval: G::ScalarField,
}
