use super::LigeroPCParams;
use super::LinCodeInfo;
use crate::{PCCommitterKey, PCUniversalParams, PCVerifierKey};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;

impl<F, C, H> LigeroPCParams<F, C, H>
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

impl<F, C, H> PCUniversalParams for LigeroPCParams<F, C, H>
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

impl<F, C, H> PCCommitterKey for LigeroPCParams<F, C, H>
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
        <LigeroPCParams<F, C, H> as PCCommitterKey>::max_degree(self)
    }
}

impl<F, C, H> PCVerifierKey for LigeroPCParams<F, C, H>
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
        <LigeroPCParams<F, C, H> as PCVerifierKey>::max_degree(self)
    }
}

impl<F, C, H> LinCodeInfo<C, H> for LigeroPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn check_well_formedness(&self) -> bool {
        self.check_well_formedness
    }

    fn distance(&self) -> (usize, usize) {
        (self.rho_inv - 1, self.rho_inv)
    }

    fn sec_param(&self) -> usize {
        self.sec_param
    }

    fn leaf_hash_params(&self) -> &<<C as Config>::LeafHash as CRHScheme>::Parameters {
        &self.leaf_hash_params
    }

    fn two_to_one_params(&self) -> &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters {
        &self.two_to_one_params
    }

    fn col_hash_params(&self) -> &<H as CRHScheme>::Parameters {
        &self.col_hash_params
    }
}
