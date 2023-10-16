use super::BreakdownPCParams;
use super::LinCodeInfo;
use crate::utils::{ceil_mul, ent};
use crate::{PCCommitterKey, PCUniversalParams, PCVerifierKey};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use ark_ff::PrimeField;
use ark_std::vec::Vec;
#[cfg(not(feature = "std"))]
use num_traits::Float;

impl<F, C, H> PCUniversalParams for BreakdownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        todo!()
    }
}

impl<F, C, H> PCCommitterKey for BreakdownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        <BreakdownPCParams<F, C, H> as PCCommitterKey>::max_degree(self)
    }
}

impl<F, C, H> PCVerifierKey for BreakdownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        <BreakdownPCParams<F, C, H> as PCVerifierKey>::max_degree(self)
    }
}

impl<F, C, H> LinCodeInfo<C, H> for BreakdownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
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

    fn col_hash_params(&self) -> &<H as CRHScheme>::Parameters {
        &self.col_hash_params
    }
}

impl<F, C, H> BreakdownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
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
    /// cn_const
    fn cn_const(&self) -> (f64, f64) {
        let a = div(self.alpha);
        let b = div(self.beta);
        let arg = 1.28 * b / a;
        let nom = ent(b) + a * ent(arg);
        let den = -b * arg.log2();
        (nom, den)
    }
    /// cn
    fn cn(&self, n: usize) -> usize {
        use ark_std::cmp::{max, min};
        let b = self.beta;
        let c = self.cn_const();
        min(
            max(ceil_mul(n, (32 * b.0, 25 * b.1)), 4 + ceil_mul(n, b)),
            ((110f64 / (n as f64) + c.0) / c.1).ceil() as usize,
        )
    }
    /// dn_const
    fn dn_const(&self) -> (f64, f64) {
        let a = div(self.alpha);
        let b = div(self.beta);
        let r = div(self.rho_inv);
        let m = self.mu();
        let n = self.nu();
        let nm = n / m;
        let nom = r * a * ent(b / r) + m * ent(nm);
        let den = -a * b * nm.log2();
        (nom, den)
    }
    /// dn
    fn dn(&self, n: usize) -> usize {
        let b = self.beta;
        let r = self.rho_inv;
        let d = self.dn_const();
        let v1 = {
            ceil_mul(n, (2 * b.0, b.1)) + // 2 * beta * n 
            ((ceil_mul(n, r) - n + 110) // n * (r - 1 + 110/n)
            as f64 / F::MODULUS_BIT_SIZE as f64).ceil() as usize
        };
        let v2 = ((110f64 / (n as f64) + d.0) / d.1).ceil() as usize;
        if v1 < v2 {
            v1
        } else {
            v2
        }
    }
    fn mat_size(&self) -> (Vec<(usize, usize, usize)>, Vec<(usize, usize, usize)>) {
        assert!(self.n > self.base_len); // TODO move this to new function
        let mut a_dims: Vec<(usize, usize, usize)> = Vec::default();
        let a = self.alpha;
        let r = self.rho_inv;

        let mut n = self.n;
        while n > self.base_len {
            let m = ceil_mul(n, a);
            let cn = self.cn(n);
            let cn = if cn < m { cn } else { m }; // can't generate more nonzero entries than there are columns
            a_dims.push((n, m, cn));
            n = m;
        }

        let b_dims = a_dims
            .iter()
            .map(|&(an, am, _)| {
                let n = ceil_mul(am, r);
                let m = ceil_mul(an, r) - an - n;
                let dn = self.dn(n);
                let dn = if dn < m { dn } else { m }; // can't generate more nonzero entries than there are columns
                (n, m, dn)
            })
            .collect::<Vec<_>>();
        (a_dims, b_dims)
    }

    fn codeword_len(&self) -> usize {
        let (a_dims, b_dims) = self.mat_size();
        b_dims.iter().map(|(_, m, _)| m).sum::<usize>() + // Output v of the recursive encoding
        a_dims.iter().map(|(n, _, _)| n).sum::<usize>() + // Input x to the recursive encoding
        b_dims.last().unwrap().0 // Output z of the last step of recursion
    }
}

#[inline]
fn div(a: (usize, usize)) -> f64 {
    a.0 as f64 / a.1 as f64
}
