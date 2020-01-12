#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

/*
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(D B^T A^-1 B)
 * where D is a square matrix and the LDLT_factor of A is provided.
 */

TEST(MathMatrixPrimMat, trace_gen_inv_quad_form_ldlt) {
  using stan::math::matrix_d;
  using stan::math::trace_gen_inv_quad_form_ldlt;

  matrix_d D00(0, 0), B00(0, 0), B02(0, 2);
  stan::math::LDLT_factor<double, 0, 0> ldlt_A0;
  EXPECT_FLOAT_EQ(0, trace_gen_inv_quad_form_ldlt(D00, ldlt_A0, B00));
  EXPECT_THROW(trace_gen_inv_quad_form_ldlt(D00, ldlt_A0, B02),
               std::invalid_argument);

  matrix_d D(2, 2), A(4, 4), B(4, 2), gen_inv_quad_form;

  D << 1, 2, 3, 4;
  A << 9.0, 3.0, 3.0, 3.0, 3.0, 10.0, 2.0, 2.0, 3.0, 2.0, 7.0, 1.0, 3.0, 2.0,
      1.0, 112.0;
  B << 100, 10, 0, 1, -3, -3, 5, 2;

  gen_inv_quad_form = D * B.transpose() * A.inverse() * B;

  stan::math::LDLT_factor<double, -1, -1> ldlt_A;
  ldlt_A.compute(A);
  ASSERT_TRUE(ldlt_A.success());

  EXPECT_FLOAT_EQ(stan::math::trace(gen_inv_quad_form),
                  stan::math::trace_gen_inv_quad_form_ldlt(D, ldlt_A, B));
}
