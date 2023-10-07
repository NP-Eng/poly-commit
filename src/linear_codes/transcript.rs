use core::marker::PhantomData;

use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use merlin::Transcript;

use crate::{to_bytes, Error};

/// The following struct is taken from jellyfish repository. Once they change
/// their dependency on `crypto-primitive`, we use their crate instead of
/// a copy-paste. We needed the newer `crypto-primitive` for serializing.
#[derive(Clone)]
pub(crate) struct IOPTranscript<F: PrimeField> {
    transcript: Transcript,
    is_empty: bool,
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

// TODO: merge this with jf_plonk::transcript
impl<F: PrimeField> IOPTranscript<F> {
    /// Create a new IOP transcript.
    pub(crate) fn new(label: &'static [u8]) -> Self {
        Self {
            transcript: Transcript::new(label),
            is_empty: true,
            phantom: PhantomData,
        }
    }

    /// Append the message to the transcript.
    pub(crate) fn append_message(&mut self, label: &'static [u8], msg: &[u8]) -> Result<(), Error> {
        self.transcript.append_message(label, msg);
        self.is_empty = false;
        Ok(())
    }

    /// Append the message to the transcript.
    pub(crate) fn append_serializable_element<S: CanonicalSerialize>(
        &mut self,
        label: &'static [u8],
        group_elem: &S,
    ) -> Result<(), Error> {
        self.append_message(
            label,
            &to_bytes!(group_elem).map_err(|_| Error::TranscriptError)?,
        )
    }

    /// Generate the challenge from the current transcript
    /// and append it to the transcript.
    ///
    /// The output field element is statistical uniform as long
    /// as the field has a size less than 2^384.
    pub(crate) fn get_and_append_challenge(&mut self, label: &'static [u8]) -> Result<F, Error> {
        //  we need to reject when transcript is empty
        if self.is_empty {
            return Err(Error::TranscriptError);
        }

        let mut buf = [0u8; 64];
        self.transcript.challenge_bytes(label, &mut buf);
        let challenge = F::from_le_bytes_mod_order(&buf);
        self.append_serializable_element(label, &challenge)?;
        Ok(challenge)
    }

    /// Generate the challenge from the current transcript
    /// and append it to the transcript.
    ///
    /// Without exposing the internal field `transcript`,
    /// this is a wrapper around getting bytes as opposed to field elements.
    pub(crate) fn get_and_append_byte_challenge(
        &mut self,
        label: &'static [u8],
        dest: &mut [u8],
    ) -> Result<(), Error> {
        //  we need to reject when transcript is empty
        if self.is_empty {
            return Err(Error::TranscriptError);
        }

        self.transcript.challenge_bytes(label, dest);
        self.append_message(label, dest)?;
        Ok(())
    }
}
