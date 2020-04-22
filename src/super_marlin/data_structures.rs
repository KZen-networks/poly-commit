use crate::*;
use algebra_core::{AffineCurve, PairingEngine, PrimeField, ProjectiveCurve, ToBytes, Zero};
use core::ops::{Add, AddAssign};

use class_group::primitives::polynomial_comm::*;
use class_group::BinaryQF;
use curv::BigInt;
use curv::FE;
use curv::elliptic::curves::traits::ECScalar;
use curv::arithmetic::traits::Converter;
use std::ops::Div;
/// `UniversalParams` are the universal parameters for the KZG10 scheme.
#[derive(PartialEq, Eq,Clone, Debug)]
pub struct UniversalParams {
    pub d_max: usize,
    pub pp: PP,


}

impl PCUniversalParams for UniversalParams {
    fn max_degree(&self) -> usize {
        self.d_max
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct CommitterKey {
    pub uni_params : UniversalParams,
    pub supported_degree: usize,
}

impl PCCommitterKey for CommitterKey{
    fn max_degree(&self) -> usize{
        self.uni_params.d_max
    }

    fn supported_degree(&self) -> usize{
        self.supported_degree
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct VerifierKey{
    pub uni_params : UniversalParams,
    pub supported_degree: usize,

}

impl PCVerifierKey for VerifierKey{
    fn max_degree(&self) -> usize{
        self.uni_params.d_max
    }

    fn supported_degree(&self) -> usize{
        self.supported_degree

    }
}

/// `Commitment` commits to a polynomial. It is output by `KZG10::commit`.
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Commitment{
    /// The commitment is a group element.
    pub comm: PolyComm,
}
/*
impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, writer: W) -> algebra_core::io::Result<()> {
        self.0.write(writer)
    }
}
impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.comm.write(&mut writer)?;
        let shifted_exists = self.shifted_comm.is_some();
        shifted_exists.write(&mut writer)?;
        self.shifted_comm
            .as_ref()
            .unwrap_or(&kzg10::Commitment::empty())
            .write(&mut writer)
    }
}
*/

impl ToBytes for Commitment {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        let bn_a= &self.comm.c.a;
        let mut bytes_a = BigInt::to_vec(&bn_a);
        let bn_b= &self.comm.c.b;
        let bytes_b = BigInt::to_vec(&bn_b);
        let bn_c= &self.comm.c.c;
        let bytes_c = BigInt::to_vec(&bn_c);
        bytes_a.extend_from_slice(&bytes_b);
        bytes_a.extend_from_slice(&bytes_c);
        Ok(writer.write_all(&bytes_a[..])?)
    }
}

impl PCCommitment for Commitment {
    #[inline]
    fn empty() -> Self {
        let bqf_empty = BinaryQF{
            a: BigInt::zero(),
            b: BigInt::zero(),
            c: BigInt::zero(),
        };
        let poly_com = PolyComm{c:bqf_empty};
        Commitment{comm: poly_com}
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        let bit_len = self.comm.c.a.bit_length() + self.comm.c.b.bit_length() + self.comm.c.c.bit_length();
        let car = if bit_len % 8 > 0 {1} else {0};
        (bit_len.div(8) + car)
     //   algebra_core::to_bytes![E::G1Affine::zero()].unwrap().len() / 2
    }
}

/*
impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Commitment<E>)> for Commitment<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Commitment<E>)) {
        let mut other = other.0.mul(f.into_repr());
        other.add_assign_mixed(&self.0);
        self.0 = other.into();
    }
}
*/

/// `Randomness` hides the polynomial inside a commitment. It is output by `KZG10::commit`.
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Randomness {
    pub f_q: BigInt,
}
/*
impl<E: PairingEngine> Randomness<E> {
    /// Does `self` provide any hiding properties to the corresponding commitment?
    /// `self.is_hiding() == true` only if the underlying polynomial is non-zero.
    #[inline]
    pub fn is_hiding(&self) -> bool {
        !self.blinding_polynomial.is_zero()
    }

    /// What is the degree of the hiding polynomial for a given hiding bound?
    #[inline]
    pub fn calculate_hiding_polynomial_degree(hiding_bound: usize) -> usize {
        hiding_bound + 1
    }
}
*/

impl PCRandomness for Randomness {
    fn empty() -> Self {
        Self {
            f_q: BigInt::zero(),
        }
    }

    //TODO: we are not really going to use it
    fn rand<R: RngCore>(hiding_bound: usize, _: bool, rng: &mut R) -> Self {
        let randomness = Randomness::empty();
        randomness
    }
}

/*

impl<'a, E: PairingEngine> Add<&'a Randomness<E>> for Randomness<E> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: &'a Self) -> Self {
        self.blinding_polynomial += &other.blinding_polynomial;
        self
    }
}

impl<'a, E: PairingEngine> Add<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: (E::Fr, &'a Randomness<E>)) -> Self {
        self += other;
        self
    }
}

impl<'a, E: PairingEngine> AddAssign<&'a Randomness<E>> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.blinding_polynomial += &other.blinding_polynomial;
    }
}

impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Randomness<E>)) {
        self.blinding_polynomial += (f, &other.blinding_polynomial);
    }
}
*/

/// `Proof` is an evaluation proof that is output by `KZG10::open`.
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Proof {
    pub proof_vec: Vec<NiEvalProof>,
}


impl PCProof for Proof {
    fn size_in_bytes(&self) -> usize {
        0 //TODO: implement
    }
}

impl ToBytes for Proof {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, writer: W) -> algebra_core::io::Result<()> {
    Ok(())

    }
}


