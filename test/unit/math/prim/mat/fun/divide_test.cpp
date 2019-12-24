#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixPrimMat, divide) {
  using stan::math::divide;
  stan::math::vector_d v0(3);
  stan::math::row_vector_d rv0(3);
  stan::math::matrix_d m0(3, 3);

  v0 << 1.0, 2.0, 3.0;
  rv0 << 1.0, 2.0, 3.0;
  m0 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;

  stan::math::vector_d v0_over_2(3);
  stan::math::row_vector_d rv0_over_2(3);
  stan::math::matrix_d m0_over_2(3, 3);
  v0_over_2 << 0.5, 1.0, 1.5;
  rv0_over_2 << 0.5, 1.0, 1.5;
  m0_over_2 << 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5;

  expect_matrix_eq(v0_over_2, divide(v0, 2.0));
  expect_matrix_eq(rv0_over_2, divide(rv0, 2.0));
  expect_matrix_eq(m0_over_2, divide(m0, 2.0));

  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (double value : {inf, -inf, nan}) {
    EXPECT_NO_THROW(divide(v0, value));
    EXPECT_NO_THROW(divide(rv0, value));
    EXPECT_NO_THROW(divide(m0, value));
    stan::math::vector_d v1 = v0;
    stan::math::row_vector_d rv1 = rv0;
    stan::math::matrix_d m1 = m0;
    v1(1) = value;
    rv1(1) = value;
    m1(1, 1) = value;
    EXPECT_NO_THROW(divide(v1, 2.0));
    EXPECT_NO_THROW(divide(rv1, 2.0));
    EXPECT_NO_THROW(divide(m1, 2.0));
  }
}
