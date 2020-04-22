mod data_structures;
pub use data_structures::*;
use rand_core::RngCore;
use crate::Error;
use crate::super_marlin::data_structures::UniversalParams;
use crate::PolynomialCommitment;
use algebra_core::Field;
use algebra_core::fields::FpCurv;
use class_group::primitives::polynomial_comm::*;
use class_group::BinaryQF;
use curv::{BigInt, FE};
use curv::elliptic::curves::traits::ECScalar;
use crate::super_marlin::data_structures::{CommitterKey, VerifierKey, Randomness, Commitment, Proof};
use crate::QuerySet;
use crate::data_structures::LabeledPolynomial;
use crate::data_structures::LabeledCommitment;
use crate::data_structures::PCUniversalParams;

pub struct SuperMarlinPoly {}

impl PolynomialCommitment<FpCurv> for SuperMarlinPoly {
    type UniversalParams = UniversalParams;
    type CommitterKey = CommitterKey;
    type VerifierKey = VerifierKey;
    type Commitment = Commitment;
    type Randomness = Randomness;
    type Proof = Proof;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>{
        let degree_bn = BigInt::from(max_degree.clone() as i32);
        let uni_params  = UniversalParams{
            d_max: max_degree,
            pp: PolyComm::setup(&degree_bn)
        };
        Ok(uni_params)

    }


    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error>{
        // lazy check that supported degree is also the max degree - otherwise return error!
        let max_degree = pp.max_degree();
        if &supported_degree > &max_degree {
            return Err(Error::TrimmingDegreeTooLarge);
        }
        //TODO:  check enforced_degree_bounds

        let mut coef_vec: Vec<FE> = Vec::new();
        let mut i = 0;
        while i < supported_degree {
            // TODO: check that i < d_max
            coef_vec.push(FE::new_random());
            i = i + 1;
        }
        let ck = CommitterKey {
            uni_params : pp.clone(),
            supported_degree ,

        };

        let vk = VerifierKey {
            uni_params : pp.clone(),
            supported_degree ,
        };
        Ok((ck, vk))

    }


    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, FpCurv>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >{
        let mut poly_com_vec: Vec<LabeledCommitment<Commitment>> = Vec::new();
        let mut randomness_vec : Vec<Randomness> = Vec::new();
        for p in polynomials{
            let coef_vec = fp_curv_elements_into_fe_elements(&p.polynomial().coeffs);
            let (comm, rand) = PolyComm::commit(&ck.uni_params.pp, &coef_vec);
            let labeled_comm = LabeledCommitment::new(
                p.label().clone(),
                Commitment{comm},
                Some(ck.supported_degree.clone()));

            poly_com_vec.push(labeled_comm);
            let randomness = Randomness{f_q: rand};
            randomness_vec.push(randomness);

        }
        Ok((poly_com_vec, randomness_vec))



    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, FpCurv>>,
        point: FpCurv,
        opening_challenge: FpCurv,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::Proof, Self::Error>
        where
            Self::Randomness: 'a{
        let mut y_vec : Vec<FpCurv> = Vec::new();
        let mut coef_vec: Vec<Vec<FE>> = Vec::new();
        for p in labeled_polynomials{
            y_vec.push(p.evaluate(point.clone()));
            coef_vec.push(fp_curv_elements_into_fe_elements(&p.polynomial().coeffs))
        }

        let proof_vec = (0..y_vec.len()).map(|i|{
            let (comm, _) = PolyComm::commit(&ck.uni_params.pp,
                                        &coef_vec[i]);
            comm.eval_prove(&ck.uni_params.pp, &point.0, &y_vec[i].0, &coef_vec[i])
        }).collect::<Vec<NiEvalProof>>();

        Ok(Proof{
            proof_vec
        })

    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: FpCurv,
        values: impl IntoIterator<Item = FpCurv>,
        proof: &Self::Proof,
        opening_challenge: FpCurv,
    ) -> Result<bool, Self::Error>
        where
            Self::Commitment: 'a{
      //  if values.len() != commitments.len(){
      //      Err(()) //TODO: add more checks
      //  }
        let mut flag = true;

        let mut i = 0;
            let iter1 = commitments.into_iter().zip(values);
            for (comm, value) in iter1 {
                let comm_c = &comm.commitment().comm;
                let proof = proof.proof_vec[i].clone();//TODO: get rid of clone
                i = i + 1;
                let y = &value.0;
                let z = point.0;
                let res = proof.eval_verify(comm_c.c.clone(), &vk.uni_params.pp, &z, y);
                if res.is_err() {
                    flag = false
                }
        }
        match flag{
            true => Ok(true),
            false => Err(Error::DegreeIsZero ),
        }

    }
}




fn fp_curv_elements_into_fe_elements(fp_curv_vec: &Vec<FpCurv>) -> Vec<FE>{
    fp_curv_vec.into_iter().map(|elem|{
        elem.0
    }).collect::<Vec<FE>>()
}



#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::super_marlin::SuperMarlinPoly;
    use algebra_core::fields::FpCurv;


    type PC = SuperMarlinPoly;


    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<FpCurv, PC>().expect("test failed for fp_curv_single_poly");

    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<FpCurv, PC>()
            .expect("test failed for fp_curv quadratic poly");

    }

        #[test]
        fn linear_poly_degree_bound_test() {
            use crate::tests::*;
            linear_poly_degree_bound_test::<FpCurv, PC>().expect("test failed for fp_curv linear");

        }

        #[test]
        fn single_poly_degree_bound_test() {
            use crate::tests::*;
            single_poly_degree_bound_test::<FpCurv, PC>().expect("test failed for fp_curv degree bound");

        }

        #[test]
        fn single_poly_degree_bound_multiple_queries_test() {
            use crate::tests::*;
            single_poly_degree_bound_multiple_queries_test::<FpCurv, PC>()
                .expect("test failed for fp_curv single_poly_degree_bound_multiple_queries");
        }

        #[test]
        fn two_polys_degree_bound_single_query_test() {
            use crate::tests::*;
            two_polys_degree_bound_single_query_test::<FpCurv, PC>()
                .expect("test failed for fp_curv");

        }

        #[test]
        fn full_end_to_end_test() {
            use crate::tests::*;
            full_end_to_end_test::<FpCurv, PC>().expect("test failed for fp_curv ");

        }

        #[test]
        fn single_equation_test() {
            use crate::tests::*;
            single_equation_test::<FpCurv, PC>().expect("test failed for fp_curv ");

        }

        #[test]
        fn two_equation_test() {
            use crate::tests::*;
            two_equation_test::<FpCurv, PC>().expect("test failed for fp_curv ");

        }

        #[test]
        fn two_equation_degree_bound_test() {
            use crate::tests::*;
            two_equation_degree_bound_test::<FpCurv, PC>().expect("test failed for fp_curv ");

        }

        #[test]
        fn full_end_to_end_equation_test() {
            use crate::tests::*;
            full_end_to_end_equation_test::<FpCurv, PC>().expect("test failed for fp_curv ");

        }

}
