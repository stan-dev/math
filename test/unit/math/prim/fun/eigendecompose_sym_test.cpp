#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, eigendecompose_sym) {
  stan::math::matrix_d m0;
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_d ev_m1(1, 1);
  ev_m1 << 2.0;

  using stan::math::eigendecompose_sym;
  EXPECT_NO_THROW(eigendecompose_sym(m0));
  EXPECT_NO_THROW(eigendecompose_sym(ev_m1));
  EXPECT_THROW(eigendecompose_sym(m1), std::invalid_argument);
}
