use core::marker::PhantomData;

use ark_ff::BigInteger;

use ark_ff::PrimeField;

use ark_std::log2;
use halo2_base::halo2_proofs::halo2curves::group::ff::PrimeField as _;
use halo2_base::utils::ScalarField;

use halo2_base::halo2_proofs::halo2curves::bn256::Fr;

use poseidon_native::Poseidon;

use crate::Error;

#[derive(Clone)]
pub(crate) struct IOPTranscript<F>
where
    F: PrimeField,
{
    pub(crate) sponge: Poseidon<Fr, 3, 2>,
    phantom: PhantomData<F>,
}

impl<F> IOPTranscript<F>
where
    F: PrimeField,
{
    pub(crate) fn new(_label: &'static [u8]) -> Self {
        Self {
            sponge: Poseidon::<Fr, 3, 2>::new(8, 57),
            phantom: PhantomData,
        }
    }

    fn absorb(&mut self, elements: &[F]) {
        for elem in elements {
            let bytes = BigInteger::to_bytes_le(&elem.into_bigint());
            let halo2_base_elem = Fr::from_bytes_le(&bytes);
            self.sponge.update(&[halo2_base_elem]);
        }
    }

    fn squeeze(&mut self) -> F {
        let halo2_out = self.sponge.squeeze();
        let out = F::from_le_bytes_mod_order(&halo2_out.to_bytes());
        out
    }

    pub(crate) fn get_and_append_challenge(&mut self, _label: &'static [u8]) -> Result<F, Error> {
        Ok(self.squeeze())
    }

    pub(crate) fn get_and_append_byte_challenge(
        &mut self,
        _label: &'static [u8],
        n: usize,
    ) -> Result<Vec<u8>, Error> {
        let out = self.sponge.squeeze();
        let bits = out
            .to_u64_limbs(Fr::NUM_BITS as usize, 1)
            .into_iter()
            .map(|x| {
                // println!("x: {}", x);
                x as u8
            })
            .collect::<Vec<u8>>();
        let bits = bits[..(log2(n) as usize)].to_vec();
        Ok(bits)
    }

    /// Append the message to the transcript.
    pub(crate) fn append_field_elements(
        &mut self,
        _label: &'static [u8],
        group_elem: &[F],
    ) -> Result<(), Error> {
        self.absorb(group_elem);
        Ok(())
    }

    /// Append the message to the transcript.
    pub(crate) fn append_serializable_element(
        &mut self,
        _label: &'static [u8],
        group_elem: &F,
    ) -> Result<(), Error> {
        self.absorb(&[*group_elem]);
        Ok(())
    }
}
