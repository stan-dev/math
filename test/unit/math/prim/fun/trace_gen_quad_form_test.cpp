#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, trace_gen_quad_form) {
  using stan::math::matrix_d;
  using stan::math::trace_gen_quad_form;

  matrix_d ad(4, 4);
  matrix_d bd(4, 2);
  matrix_d cd(2, 2);
  double res;

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;
  cd.setIdentity(2, 2);

  res = trace_gen_quad_form(cd, ad, bd);
  EXPECT_FLOAT_EQ(26758, res);

  EXPECT_THROW(trace_gen_quad_form(ad, ad, bd), std::invalid_argument);
  EXPECT_THROW(trace_gen_quad_form(bd, ad, bd), std::invalid_argument);
  EXPECT_THROW(trace_gen_quad_form(ad, bd, bd), std::invalid_argument);
}

TEST(MathMatrixPrim, trace_gen_quad_form_size_zero) {
  using stan::math::matrix_d;
  using stan::math::trace_gen_quad_form;

  matrix_d a00, b00, d00;
  matrix_d b02(0, 2);
  matrix_d d22(2, 2);
  d22 << 1, 2, 3, 4;
  double res;

  res = trace_gen_quad_form(d00, a00, b00);
  EXPECT_FLOAT_EQ(0, res);

  res = trace_gen_quad_form(d22, a00, b02);
  EXPECT_FLOAT_EQ(0, res);
}
