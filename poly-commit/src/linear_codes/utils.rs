use crate::{utils::ceil_div, Error};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::string::ToString;
use ark_std::{collections::BTreeSet, vec::Vec};

#[cfg(all(not(feature = "std"), target_arch = "aarch64"))]
use num_traits::Float;

#[cfg(test)]
use {
    ark_crypto_primitives::crh::CRHScheme,
    ark_std::{borrow::Borrow, rand::RngCore},
};

/// Apply reed-solomon encoding to msg.
/// Assumes msg.len() is equal to the order of some FFT domain in F.
/// Returns a vector of length equal to the smallest FFT domain of size at least msg.len() * RHO_INV.
pub(crate) fn reed_solomon<F: FftField>(
    // msg, of length m, is interpreted as a vector of coefficients of a polynomial of degree m - 1
    msg: &[F],
    rho_inv: usize,
) -> Vec<F> {
    let m = msg.len();

    let extended_domain = GeneralEvaluationDomain::<F>::new(m * rho_inv).unwrap_or_else(|| {
        panic!(
            "The field F cannot accomodate FFT for msg.len() * RHO_INV = {} elements (too many)",
            m * rho_inv
        )
    });

    extended_domain.fft(msg)
}

/// This is CSC format
/// https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct SprsMat<F: Field> {
    /// Number of rows.
    pub(crate) n: usize,
    /// Number of columns.
    pub(crate) m: usize,
    /// Number of non-zero entries in each row.
    pub(crate) d: usize,
    /// Numbers of non-zero elements in each columns.
    ind_ptr: Vec<usize>,
    /// The indices in each columns where exists a non-zero element.
    col_ind: Vec<usize>,
    // The values of non-zero entries.
    val: Vec<F>,
}

impl<F: Field> SprsMat<F> {
    /// Calulates v.M
    pub(crate) fn row_mul(&self, v: &[F]) -> Vec<F> {
        (0..self.m)
            .map(|j| {
                let ij = self.ind_ptr[j]..self.ind_ptr[j + 1];
                self.col_ind[ij.clone()]
                    .iter()
                    .zip(&self.val[ij])
                    .map(|(&idx, x)| v[idx] * x)
                    .sum::<F>()
            })
            .collect::<Vec<_>>()
    }
    /// Create a new `SprsMat` from list of elements that represents the
    /// matrix in column major order. `n` is the number of rows, `m` is
    /// the number of columns, and `d` is NNZ in each row.
    pub fn new_from_flat(n: usize, m: usize, d: usize, list: &[F]) -> Self {
        let nnz = d * n;
        let mut ind_ptr = vec![0; m + 1];
        let mut col_ind = Vec::<usize>::with_capacity(nnz);
        let mut val = Vec::<F>::with_capacity(nnz);
        assert!(list.len() == m * n, "The dimension is incorrect.");
        for i in 0..m {
            for (c, &v) in list[i * n..(i + 1) * n].iter().enumerate() {
                if v != F::zero() {
                    ind_ptr[i + 1] += 1;
                    col_ind.push(c);
                    val.push(v);
                }
            }
            ind_ptr[i + 1] += ind_ptr[i];
        }
        assert!(ind_ptr[m] <= nnz, "The dimension or NNZ is incorrect.");
        Self {
            n,
            m,
            d,
            ind_ptr,
            col_ind,
            val,
        }
    }
    pub fn new_from_columns(n: usize, m: usize, d: usize, list: &[Vec<(usize, F)>]) -> Self {
        let nnz = d * n;
        let mut ind_ptr = vec![0; m + 1];
        let mut col_ind = Vec::<usize>::with_capacity(nnz);
        let mut val = Vec::<F>::with_capacity(nnz);
        assert!(list.len() == m, "The dimension is incorrect.");
        for j in 0..m {
            for (i, v) in list[j].iter() {
                ind_ptr[j + 1] += 1;
                col_ind.push(*i);
                val.push(*v);
            }
            assert!(list[j].len() <= n, "The dimension is incorrect.");
            ind_ptr[j + 1] += ind_ptr[j];
        }
        assert!(ind_ptr[m] <= nnz, "The dimension or NNZ is incorrect.");
        Self {
            n,
            m,
            d,
            ind_ptr,
            col_ind,
            val,
        }
    }
}

#[inline]
pub(crate) fn get_num_bytes(n: usize) -> usize {
    ceil_div((usize::BITS - n.leading_zeros()) as usize, 8)
}

/// Generate `t` (not necessarily distinct) random points in `[0, n)`
/// using the current state of the `transcript`.
pub(crate) fn get_indices_from_sponge<S: CryptographicSponge>(
    n: usize,
    t: usize,
    sponge: &mut S,
) -> Result<Vec<usize>, Error> {
    let bytes_per_index = get_num_bytes(n);
    let indices_per_squeeze = 32 / bytes_per_index;
    let num_squeezes = (t + indices_per_squeeze - 1) / indices_per_squeeze;
    let mut indices = Vec::with_capacity(t);
    for i in 0..num_squeezes {
        let i_as_bytes: [u8; 2] = (i as u16).to_be_bytes();
        let tag = "index_".to_string().as_bytes().to_vec();
        let tag = [tag, i_as_bytes.to_vec()].concat();

        sponge.absorb(&tag);

        let bytes = sponge.squeeze_bytes(32);

        let bytes_to_take = if i == num_squeezes - 1 {
            t % indices_per_squeeze
        } else {
            indices_per_squeeze
        };

        let squeeze_indices = bytes
            .chunks(bytes_per_index)
            .take(bytes_to_take)
            .filter(|x| x.len() == bytes_per_index)
            .map(|x| x.iter().fold(0, |acc, &x| ((acc << 8) + x as usize) % n));

        indices.extend(squeeze_indices);
    }
    let mut set = BTreeSet::new();
    let indices = indices
        .into_iter()
        .filter(|x| set.insert(*x))
        .take(t)
        .collect::<Vec<usize>>();
    Ok(indices)
}

#[inline]
pub(crate) fn calculate_t<F: PrimeField>(
    sec_param: usize,
    distance: (usize, usize),
    codeword_len: usize,
) -> Result<usize, Error> {
    // Took from the analysis by BCI+20 and Ligero
    // We will find the smallest $t$ such that
    // $(1-\delta)^t + (\rho+\delta)^t + \frac{n}{F} < 2^{-\lambda}$.
    // With $\delta = \frac{1-\rho}{2}$, the expreesion is
    // $2 * (\frac{1+\rho}{2})^t + \frac{n}{F} < 2^(-\lambda)$.

    let field_bits = F::MODULUS_BIT_SIZE as i32;
    let sec_param = sec_param as i32;

    let residual = codeword_len as f64 / 2.0_f64.powi(field_bits);
    let rhs = (2.0_f64.powi(-sec_param) - residual).log2();
    if !(rhs.is_normal()) {
        return Err(Error::InvalidParameters("For the given codeword length and the required security guarantee, the field is not big enough.".to_string()));
    }
    let nom = rhs - 1.0;
    let denom = (1.0 - 0.5 * distance.0 as f64 / distance.1 as f64).log2();
    if !(denom.is_normal()) {
        return Err(Error::InvalidParameters(
            "The distance is wrong".to_string(),
        ));
    }
    let t = (nom / denom).ceil() as usize;
    Ok(if t < codeword_len { t } else { codeword_len })
}

#[cfg(test)]
pub(crate) struct LeafIdentityHasher;

#[cfg(test)]
impl CRHScheme for LeafIdentityHasher {
    type Input = Vec<u8>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        Ok(input.borrow().to_vec().into())
    }
}

pub(crate) fn tensor_vec<F: PrimeField>(values: &[F]) -> Vec<F> {
    let one = F::one();
    let anti_values: Vec<F> = values.iter().map(|v| one - *v).collect();

    let mut layer: Vec<F> = vec![one];

    for i in 0..values.len() {
        let mut new_layer = Vec::new();
        for v in &layer {
            new_layer.push(*v * anti_values[i]);
        }
        for v in &layer {
            new_layer.push(*v * values[i]);
        }
        layer = new_layer;
    }

    layer
}

#[cfg(test)]
pub(crate) mod tests {

    use crate::utils::to_field;

    use super::*;

    use ark_bls12_377::Fq;
    use ark_bls12_377::Fr;
    use ark_crypto_primitives::sponge::Absorb;
    use ark_crypto_primitives::sponge::FieldElementSize;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        Polynomial,
    };
    use ark_std::test_rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
    use tiny_keccak::Hasher;
    use tiny_keccak::Keccak;

    #[derive(Clone)]
    struct TestSponge {
        pub(crate) state: Vec<u8>,
    }

    trait CustomCryptographicSponge: CryptographicSponge {
        fn set_state(&mut self, state: Vec<u8>);
    }

    impl CryptographicSponge for TestSponge {
        type Config = ();

        fn new(_params: &Self::Config) -> Self {
            TestSponge { state: vec![] }
        }

        fn absorb(&mut self, input: &impl Absorb) {
            let mut input_bytes = vec![];
            input.to_sponge_bytes(&mut input_bytes);
            self.state.extend_from_slice(&input_bytes);
        }

        fn squeeze_bytes(&mut self, num_bytes: usize) -> Vec<u8> {
            let mut keccak = Keccak::v256();
            let mut output = vec![0u8; num_bytes];
            keccak.update(&self.state);
            keccak.finalize(&mut output);
            self.state = output.clone();
            output
        }

        fn squeeze_bits(&mut self, _num_bits: usize) -> Vec<bool> {
            unimplemented!("squeeze_bits is not implemented for TestSponge")
        }

        fn squeeze_field_elements_with_sizes<F: PrimeField>(
            &mut self,
            sizes: &[FieldElementSize],
        ) -> Vec<F> {
            unimplemented!("squeeze_field_elements_with_sizes is not implemented for TestSponge")
        }
    }

    fn test_sponge() -> TestSponge {
        TestSponge::new(&())
    }

    #[test]
    fn test_reed_solomon() {
        let rho_inv = 3;
        // `i` is the min number of evaluations we need to interpolate a poly of degree `i - 1`
        for i in 1..10 {
            let deg = (1 << i) - 1;

            let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
            let mut pol = DensePolynomial::rand(deg, rand_chacha);

            while pol.degree() != deg {
                pol = DensePolynomial::rand(deg, rand_chacha);
            }

            let coeffs = &pol.coeffs;

            // size of evals might be larger than deg + 1 (the min. number of evals needed to interpolate): we could still do R-S encoding on smaller evals, but the resulting polynomial will differ, so for this test to work we should pass it in full
            let m = deg + 1;

            let encoded = reed_solomon(&coeffs, rho_inv);

            let large_domain = GeneralEvaluationDomain::<Fr>::new(m * rho_inv).unwrap();

            // the encoded elements should agree with the evaluations of the polynomial in the larger domain
            for j in 0..(rho_inv * m) {
                assert_eq!(pol.evaluate(&large_domain.element(j)), encoded[j]);
            }
        }
    }

    #[test]
    fn test_sprs_row_mul() {
        // The columns major representation of a matrix.
        let mat: Vec<Fr> = to_field(vec![10, 23, 55, 100, 1, 58, 4, 0, 9]);

        let mat = SprsMat::new_from_flat(3, 3, 3, &mat);
        let v: Vec<Fr> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fr
        // we ensure the test will still pass even if Fr changes
        assert_eq!(mat.row_mul(&v), to_field::<Fr>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_sprs_row_mul_sparse_mat() {
        // The columns major representation of a matrix.
        let mat: Vec<Fr> = to_field(vec![10, 23, 55, 100, 1, 58, 4, 0, 9]);
        let mat = vec![
            vec![(0usize, mat[0]), (1usize, mat[1]), (2usize, mat[2])],
            vec![(0usize, mat[3]), (1usize, mat[4]), (2usize, mat[5])],
            vec![(0usize, mat[6]), (1usize, mat[7]), (2usize, mat[8])],
        ];

        let mat = SprsMat::new_from_columns(3, 3, 3, &mat);
        let v: Vec<Fr> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fr
        // we ensure the test will still pass even if Fr changes
        assert_eq!(mat.row_mul(&v), to_field::<Fr>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_get_num_bytes() {
        assert_eq!(get_num_bytes(0), 0);
        assert_eq!(get_num_bytes(1), 1);
        assert_eq!(get_num_bytes(9), 1);
        assert_eq!(get_num_bytes(1 << 11), 2);
        assert_eq!(get_num_bytes(1 << 32 - 1), 4);
        assert_eq!(get_num_bytes(1 << 32), 5);
        assert_eq!(get_num_bytes(1 << 32 + 1), 5);
    }

    #[test]
    fn test_calculate_t_with_good_parameters() {
        assert!(calculate_t::<Fq>(128, (3, 4), 2_usize.pow(32)).unwrap() < 200);
        assert!(calculate_t::<Fq>(256, (3, 4), 2_usize.pow(32)).unwrap() < 400);
    }

    #[test]
    fn test_calculate_t_with_bad_parameters() {
        calculate_t::<Fq>(
            (Fq::MODULUS_BIT_SIZE - 60) as usize,
            (3, 4),
            2_usize.pow(60),
        )
        .unwrap_err();
        calculate_t::<Fq>(400, (3, 4), 2_usize.pow(32)).unwrap_err();
    }

    #[test]
    fn test_get_indices_from_sponge() {
        let expected_indices: Vec<usize> = vec![
            12828, 10294, 5381, 4213, 2882, 4840, 16057, 6998, 649, 485, 5488, 10443, 9616, 6686,
            8535, 1221, 6420, 11745, 5346, 11735, 9150, 127, 4286, 14291, 4079, 12212, 15017, 2227,
            13677, 9465, 993, 5775, 8407, 3513, 5573, 15504, 834, 8782, 12879, 6655, 2583, 3490,
            589, 5376, 5677, 12096, 12047, 2821, 15565, 6221, 10275, 1528, 12274, 819, 14782, 6792,
            116, 3241, 5430, 4516, 3339, 935, 4125, 7446, 999, 14910, 5166, 9430, 11872, 9944,
            3104, 4597, 14666, 57, 7824, 1599, 12663, 2079, 11938, 10533, 13653, 12674, 7435, 4997,
            6673, 10856, 13988, 5413, 14721, 8174, 12869, 13075, 12398, 7079, 3672, 10020, 12003,
            2988, 7038, 6553, 9777, 9533, 7171, 9530, 11512, 16147, 9769, 8116, 3703, 2758, 9342,
            2382, 13165, 11855, 12514, 4396, 910, 15236, 1079, 4606, 12979, 9489, 1310, 343, 1930,
            8772, 3418, 13781, 3541, 3485, 8599, 15356, 15457, 13185, 10404, 7389, 535, 5974, 5866,
            6132, 9321, 3586, 7027, 12394, 6097, 15669, 2811, 11237, 5221, 14039, 15331, 12991,
            10820, 2638, 6677, 446, 9666, 13817, 3208, 4196, 1440, 1497, 13098, 4861, 9306, 8135,
            14593, 1272, 7798, 12619, 9930, 5663, 14015, 2485, 13160, 2534, 15355, 14777, 15539,
            5081, 8804, 1199, 12574, 12789, 9701, 15640, 2278, 7275, 115, 11158, 7382, 712, 11337,
            14868, 2576, 12443, 3353, 8358, 16337, 3008, 10849, 578, 3615, 3265, 7557, 3345, 6186,
            6267, 2789, 6094, 16197, 12200, 13326, 7272, 3700, 9594, 16341, 11324, 798, 11222,
            5390, 987, 13510, 13606, 498, 7586, 15550, 14803, 918, 1154, 1436, 6864, 10938, 4025,
            6707, 9708, 3315, 10495, 2226, 4164, 15231, 6272, 16374, 470, 11139, 592, 2490, 5447,
            12739, 9509, 9256, 11578, 1970, 90, 12346, 1881, 2540, 13989, 7663, 15811, 7532, 7860,
            15573, 12884, 2604, 6255, 15991, 7342, 3912, 14151, 10575, 374, 9997, 3895, 12253,
            12547, 13724, 14325, 11730, 10397, 15268, 6630, 2028, 6422, 6269, 16156, 4123, 15443,
            15492, 4001, 12452, 2062, 1657, 5411, 15597, 10761, 12235, 12595, 11059, 8638, 13444,
            15038, 14506, 9839, 13217, 12752, 11534, 15461, 13322,
        ];

        let mut sponge = test_sponge();
        let n_ext_cols = 16384;
        let t = 311;
        let indices = get_indices_from_sponge(n_ext_cols, t, &mut sponge).unwrap();
        assert_eq!(indices, expected_indices);
    }
}
