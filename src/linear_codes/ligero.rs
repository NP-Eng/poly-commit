use crate::{PCCommitterKey, PCUniversalParams, PCVerifierKey};
use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;

use super::LigeroPCParams;
use super::LinCodeInfo;

impl<F, C> LigeroPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    /// Create new UniversalParams
    pub fn new(
        sec_param: usize,
        rho_inv: usize,
        check_well_formedness: bool,
        leaf_hash_params: LeafParam<C>,
        two_to_one_params: TwoToOneParam<C>,
    ) -> Self {
        Self {
            _field: PhantomData,
            sec_param,
            rho_inv,
            check_well_formedness,
            leaf_hash_params,
            two_to_one_params,
        }
    }
}

impl<F, C> PCUniversalParams for LigeroPCParams<F, C>
where
    F: PrimeField,
    C: Config,
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

impl<F, C> PCCommitterKey for LigeroPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        <LigeroPCParams<F, C> as PCCommitterKey>::max_degree(self)
    }
}

impl<F, C> PCVerifierKey for LigeroPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        <LigeroPCParams<F, C> as PCVerifierKey>::max_degree(self)
    }
}

impl<F, C> LinCodeInfo<C> for LigeroPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn check_well_formedness(&self) -> bool {
        self.check_well_formedness
    }

    fn rho_inv(&self) -> (usize, usize) {
        (self.rho_inv, 1)
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
}
