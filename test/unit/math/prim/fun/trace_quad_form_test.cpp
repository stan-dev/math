#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, trace_quad_form) {
  using stan::math::matrix_d;
  using stan::math::trace_quad_form;

  matrix_d ad(4, 4);
  matrix_d bd(4, 2);
  double res;
  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;

  res = trace_quad_form(ad, bd);
  EXPECT_FLOAT_EQ(26758, res);

  EXPECT_THROW(trace_quad_form(bd, ad), std::invalid_argument);
  EXPECT_THROW(trace_quad_form(bd, bd), std::invalid_argument);
}

TEST(MathMatrixPrim, trace_quad_form_size_zero) {
  using stan::math::matrix_d;
  using stan::math::trace_quad_form;

  matrix_d a00, b00;
  matrix_d b02(0, 2);
  double res;

  res = trace_quad_form(a00, b00);
  EXPECT_FLOAT_EQ(0, res);

  res = trace_quad_form(a00, b02);
  EXPECT_FLOAT_EQ(0, res);
}
