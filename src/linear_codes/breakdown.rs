use crate::{PCCommitterKey, PCUniversalParams, PCVerifierKey};
use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use ark_ff::PrimeField;
// use ark_std::marker::PhantomData;

use super::BreakdownPCParams;
use super::LinCodeInfo;

impl<F, C> PCUniversalParams for BreakdownPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn max_degree(&self) -> usize {
        todo!()
    }
}

impl<F, C> PCCommitterKey for BreakdownPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        <BreakdownPCParams<F, C> as PCCommitterKey>::max_degree(self)
    }
}

impl<F, C> PCVerifierKey for BreakdownPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        <BreakdownPCParams<F, C> as PCVerifierKey>::max_degree(self)
    }
}

impl<F, C> LinCodeInfo<C> for BreakdownPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    fn check_well_formedness(&self) -> bool {
        self.check_well_formedness
    }

    fn rho_inv(&self) -> (usize, usize) {
        (self.rho_inv.0, self.rho_inv.1)
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

impl<F, C> BreakdownPCParams<F, C>
where
    F: PrimeField,
    C: Config,
{
    /// Create new UniversalParams
    pub fn new(
        _sec_param: usize,
        _alpha: (usize, usize),
        _beta: (usize, usize),
        _rho_inv: (usize, usize),
        _check_well_formedness: bool,
        _leaf_hash_params: LeafParam<C>,
        _two_to_one_params: TwoToOneParam<C>,
    ) -> Self {
        todo!()
        // Self {
        //     _field: PhantomData,
        //     sec_param,
        //     alpha,
        //     beta,
        //     rho_inv,
        //     check_well_formedness,
        //     leaf_hash_params,
        //     two_to_one_params,
        // }
    }
    /// mu = rho_inv - 1 - rho_inv * alpha
    fn mu(&self) -> f64 {
        let r = self.rho_inv;
        let a = self.alpha;
        let nom = r.0 * (a.1 - a.0) - r.1 * a.1;
        let den = r.1 * a.1;
        nom as f64 / den as f64
    }
    /// nu = beta + alpha * beta + 0.03
    fn nu(&self) -> f64 {
        let a = self.alpha;
        let b = self.beta;
        let c = (3usize, 100usize);
        let nom = b.0 * (a.1 + a.0) * c.1 + c.0 * b.1 * a.1;
        let den = b.1 * a.1 * c.1;
        nom as f64 / den as f64
    }
    /// Entropy function
    fn ent(x: f64) -> f64 {
        assert!(0f64 <= x && x <= 1f64);
        if x == 0f64 || x == 1f64 {
            0f64
        } else {
            -x * x.log2() - (1.0 - x) * (1.0 - x).log2()
        }
    }
}
