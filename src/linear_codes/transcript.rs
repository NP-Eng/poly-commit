use core::marker::PhantomData;
use std::fmt::{Display, Formatter, Result as FmtResult};

use ark_ff::BigInteger;

use ark_ff::PrimeField;

use ark_std::log2;
use halo2_base::halo2_proofs::halo2curves::group::ff::PrimeField as _;
use halo2_base::utils::ScalarField;

use halo2_base::halo2_proofs::halo2curves::bn256::Fr;

use poseidon_native::Poseidon;

use crate::Error;

#[derive(Clone)]
pub(crate) enum Operation {
    Update(Vec<Fr>, String),
    Squeeze(Vec<Fr>, String),
}

#[derive(Clone)]
pub(crate) struct IOPTranscript<F>
where
    F: PrimeField,
{
    pub(crate) sponge: Poseidon<Fr, 3, 2>,
    operations: Vec<Operation>,
    phantom: PhantomData<F>,
}

impl Display for Operation {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let (op, elems, tag) = match self {
            Operation::Update(elems, tag) => ("Update", elems, tag),
            Operation::Squeeze(elems, tag) => ("Squeeze", elems, tag),
        };

        writeln!(
            f,
            "{op} ({} element{}), tag: \"{tag}\"",
            elems.len(),
            if elems.len() == 1 { "" } else { "s" }
        )?;

        for e in elems {
            writeln!(f, "    {e:?}")?;
        }

        FmtResult::Ok(())
    }
}

impl<F> IOPTranscript<F>
where
    F: PrimeField,
{
    pub(crate) fn new(_label: &'static [u8]) -> Self {
        Self {
            sponge: Poseidon::<Fr, 3, 2>::new(8, 57),
            operations: Vec::new(),
            phantom: PhantomData,
        }
    }

    fn absorb(&mut self, label: &'static [u8], elements: &[F]) {
        let mut temp_halo_values = Vec::new();
        for elem in elements {
            let bytes = BigInteger::to_bytes_le(&elem.into_bigint());
            let halo2_base_elem = Fr::from_bytes_le(&bytes);
            temp_halo_values.push(halo2_base_elem);
            self.sponge.update(&[halo2_base_elem]);
        }
        self.operations.push(Operation::Update(
            temp_halo_values,
            String::from_utf8(label.to_vec()).unwrap(),
        ));
    }

    fn squeeze(&mut self) -> F {
        let halo2_out = self.sponge.squeeze();
        let out = F::from_le_bytes_mod_order(&halo2_out.to_bytes());
        self.operations
            .push(Operation::Squeeze(vec![halo2_out], "squeeze".into()));
        out
    }

    pub(crate) fn get_and_append_challenge(&mut self, _label: &'static [u8]) -> Result<F, Error> {
        Ok(self.squeeze())
    }

    pub(crate) fn get_and_append_byte_challenge(
        &mut self,
        _label: &'static [u8],
        n: usize,
        t: usize,
    ) -> Result<Vec<Vec<u8>>, Error> {
        let mut tmp_outs = Vec::new();

        let outer_bits = (0..t)
            .map(|_| {
                let out = self.sponge.squeeze();
                tmp_outs.push(out);
                let bits = out
                    .to_u64_limbs(Fr::NUM_BITS as usize, 1)
                    .into_iter()
                    .map(|x| x as u8)
                    .collect::<Vec<u8>>();
                bits[..(log2(n) as usize)].to_vec()
            })
            .collect();
        self.operations
            .push(Operation::Squeeze(tmp_outs, "column indices".into()));
        Ok(outer_bits)
    }

    /// Append the message to the transcript.
    pub(crate) fn append_field_elements(
        &mut self,
        label: &'static [u8],
        group_elem: &[F],
    ) -> Result<(), Error> {
        self.absorb(label, group_elem);
        Ok(())
    }

    /// Append the message to the transcript.
    pub(crate) fn append_serializable_element(
        &mut self,
        label: &'static [u8],
        group_elem: &F,
    ) -> Result<(), Error> {
        self.absorb(label, &[*group_elem]);
        Ok(())
    }
}

impl<F: PrimeField> Display for IOPTranscript<F> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        writeln!(f, "{} operations:", self.operations.len())?;

        for (i, e) in self.operations.iter().enumerate() {
            write!(f, " {}. {}", i + 1, e)?;
        }

        FmtResult::Ok(())
    }
}
