use super::utils::SprsMat;
use super::BrakedownPCParams;
use super::LinCodeInfo;
use crate::linear_codes::utils::calculate_t;
use crate::utils::ceil_div;
use crate::utils::{ceil_mul, ent};
use crate::{PCCommitterKey, PCUniversalParams, PCVerifierKey};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use ark_ff::PrimeField;
use ark_std::log2;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
#[cfg(not(feature = "std"))]
use num_traits::Float;

impl<F, C, H> PCUniversalParams for BrakedownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        usize::MAX
    }
}

impl<F, C, H> PCCommitterKey for BrakedownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        usize::MAX
    }

    fn supported_degree(&self) -> usize {
        <BrakedownPCParams<F, C, H> as PCCommitterKey>::max_degree(self)
    }
}

impl<F, C, H> PCVerifierKey for BrakedownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        usize::MAX
    }

    fn supported_degree(&self) -> usize {
        <BrakedownPCParams<F, C, H> as PCVerifierKey>::max_degree(self)
    }
}

impl<F, C, H> LinCodeInfo<C, H> for BrakedownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn check_well_formedness(&self) -> bool {
        self.check_well_formedness
    }

    fn distance(&self) -> (usize, usize) {
        (self.rho_inv.1 * self.beta.0, self.rho_inv.0 * self.beta.1)
    }

    fn sec_param(&self) -> usize {
        self.sec_param
    }

    fn compute_dimensions(&self, _n: usize) -> (usize, usize) {
        (self.n, self.m)
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

impl<F, C, H> BrakedownPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    /// Create a default UniversalParams, with the values from Fig. 2 from the paper.
    pub fn default<R: RngCore>(
        rng: &mut R,
        poly_len: usize,
        check_well_formedness: bool,
        leaf_hash_params: LeafParam<C>,
        two_to_one_params: TwoToOneParam<C>,
        col_hash_params: H::Parameters,
    ) -> Self {
        let sec_param = 128;
        let a = (178, 1000);
        let b = (61, 1000);
        let r = (1521, 1000);
        let base_len = 30;
        let t = calculate_t::<F>(sec_param, (b.0 * r.1, b.1 * r.0), poly_len).unwrap(); // we want to get a rough idea what t is
        let n = 1 << log2((ceil_div(2 * poly_len, t) as f64).sqrt().ceil() as usize);
        let m = ceil_div(poly_len, n);
        let c = Self::cn_const(a, b);
        let d = Self::dn_const(a, b, r);
        let ct = Constants { a, b, r, c, d };
        let (a_dims, b_dims) = Self::mat_size(m, base_len, &ct);
        let a_mats = Self::make_all(rng, &a_dims);
        let b_mats = Self::make_all(rng, &b_dims);

        Self::new(
            sec_param,
            a,
            b,
            r,
            base_len,
            n,
            m,
            a_dims,
            b_dims,
            a_mats,
            b_mats,
            check_well_formedness,
            leaf_hash_params,
            two_to_one_params,
            col_hash_params,
        )
    }

    /// This function creates a UniversalParams.
    // TODO There should be a sanity check on the input.
    pub fn new(
        sec_param: usize,
        a: (usize, usize),
        b: (usize, usize),
        r: (usize, usize),
        base_len: usize,
        n: usize,
        m: usize,
        a_dims: Vec<(usize, usize, usize)>,
        b_dims: Vec<(usize, usize, usize)>,
        a_mats: Vec<SprsMat<F>>,
        b_mats: Vec<SprsMat<F>>,
        check_well_formedness: bool,
        leaf_hash_params: LeafParam<C>,
        two_to_one_params: TwoToOneParam<C>,
        col_hash_params: H::Parameters,
    ) -> Self {
        let m_ext = if a_dims.is_empty() {
            ceil_mul(m, r)
        } else {
            Self::codeword_len(&a_dims, &b_dims)
        };
        let start = a_dims
            .iter()
            .scan(0, |acc, &(row, _, _)| {
                *acc += row;
                Some(*acc)
            })
            .collect::<Vec<_>>();
        let end = b_dims
            .iter()
            .scan(m_ext, |acc, &(_, col, _)| {
                *acc -= col;
                Some(*acc)
            })
            .collect::<Vec<_>>();

        Self {
            sec_param,
            alpha: a,
            beta: b,
            rho_inv: r,
            base_len,
            n,
            m,
            m_ext,
            a_dims,
            b_dims,
            start,
            end,
            a_mats,
            b_mats,
            check_well_formedness,
            leaf_hash_params,
            two_to_one_params,
            col_hash_params,
        }
    }
    /// mu = rho_inv - 1 - rho_inv * alpha
    fn mu(a: (usize, usize), r: (usize, usize)) -> f64 {
        let nom = r.0 * (a.1 - a.0) - r.1 * a.1;
        let den = r.1 * a.1;
        nom as f64 / den as f64
    }
    /// nu = beta + alpha * beta + 0.03
    fn nu(a: (usize, usize), b: (usize, usize)) -> f64 {
        let c = (3usize, 100usize);
        let nom = b.0 * (a.1 + a.0) * c.1 + c.0 * b.1 * a.1;
        let den = b.1 * a.1 * c.1;
        nom as f64 / den as f64
    }
    /// cn_const
    fn cn_const(a: (usize, usize), b: (usize, usize)) -> (f64, f64) {
        let a = div(a);
        let b = div(b);
        let arg = 1.28 * b / a;
        let nom = ent(b) + a * ent(arg);
        let den = -b * arg.log2();
        (nom, den)
    }
    /// cn
    fn cn(n: usize, ct: &Constants) -> usize {
        use ark_std::cmp::{max, min};
        let b = ct.b;
        let c = ct.c;
        min(
            max(ceil_mul(n, (32 * b.0, 25 * b.1)), 4 + ceil_mul(n, b)),
            ((110f64 / (n as f64) + c.0) / c.1).ceil() as usize,
        )
    }
    /// dn_const
    fn dn_const(a: (usize, usize), b: (usize, usize), r: (usize, usize)) -> (f64, f64) {
        let m = Self::mu(a, r);
        let n = Self::nu(a, b);
        let a = div(a);
        let b = div(b);
        let r = div(r);
        let nm = n / m;
        let nom = r * a * ent(b / r) + m * ent(nm);
        let den = -a * b * nm.log2();
        (nom, den)
    }
    /// dn
    fn dn(n: usize, ct: &Constants) -> usize {
        let b = ct.b;
        let r = ct.r;
        let d = ct.d;
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
    fn mat_size(
        mut n: usize,
        base_len: usize,
        ct: &Constants,
    ) -> (Vec<(usize, usize, usize)>, Vec<(usize, usize, usize)>) {
        let mut a_dims: Vec<(usize, usize, usize)> = Vec::default();
        let a = ct.a;
        let r = ct.r;

        while n >= base_len {
            let m = ceil_mul(n, a);
            let cn = Self::cn(n, ct);
            let cn = if cn < m { cn } else { m }; // can't generate more nonzero entries than there are columns
            a_dims.push((n, m, cn));
            n = m;
        }

        let b_dims = a_dims
            .iter()
            .map(|&(an, am, _)| {
                let n = ceil_mul(am, r);
                let m = ceil_mul(an, r) - an - n;
                let dn = Self::dn(n, ct);
                let dn = if dn < m { dn } else { m }; // can't generate more nonzero entries than there are columns
                (n, m, dn)
            })
            .collect::<Vec<_>>();
        (a_dims, b_dims)
    }

    /// This function computes the codeword length
    /// Notice that it assumes the input is bigger than base_len (i.e., a_dim is not empty)
    pub(crate) fn codeword_len(
        a_dims: &[(usize, usize, usize)],
        b_dims: &[(usize, usize, usize)],
    ) -> usize {
        b_dims.iter().map(|(_, m, _)| m).sum::<usize>() + // Output v of the recursive encoding
        a_dims.iter().map(|(n, _, _)| n).sum::<usize>() + // Input x to the recursive encoding
        b_dims.last().unwrap().0 // Output z of the last step of recursion
    }

    fn make_mat<R: RngCore>(n: usize, m: usize, d: usize, rng: &mut R) -> SprsMat<F> {
        let mut array: Vec<usize> = (0..m).collect();
        let mut mat = vec![F::zero(); n * m]; // TODO is this ok?
        for i in 0..n {
            let idxs = {
                (0..d)
                    .map(|j| {
                        let r = rng.next_u64() as usize % (m - j);
                        array.swap(r, m - 1 - j);
                        array[m - 1 - j]
                    })
                    .collect::<Vec<usize>>()
            };
            for j in idxs {
                // This puts columns together
                mat[i + n * j] = loop {
                    let r = F::rand(rng);
                    if r != F::zero() {
                        break r; // Break out of the loop and assign the non-zero value to mat[i + n * j]
                    }
                };
            }
        }
        // Notice that it is transposed now
        SprsMat::<F>::new_from_flat(n, m, d, &mat)
    }

    fn make_all<R: RngCore>(rng: &mut R, dims: &[(usize, usize, usize)]) -> Vec<SprsMat<F>> {
        dims.iter()
            .map(|(n, m, d)| Self::make_mat(*n, *m, *d, rng))
            .collect::<Vec<_>>()
    }
}

#[inline]
fn div(a: (usize, usize)) -> f64 {
    a.0 as f64 / a.1 as f64
}

struct Constants {
    a: (usize, usize),
    b: (usize, usize),
    r: (usize, usize),
    c: (f64, f64),
    d: (f64, f64),
}
