#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, matrix_exp_2x2_2x2_1) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(2, 2), m2(2, 2);

  m1 << 3, 0, 0, 4;
  m2 << exp(3), 0, 0, exp(4);

  expect_matrix_eq(m2, stan::math::matrix_exp_2x2(m1));
}

TEST(MathMatrixPrimMat, matrix_exp_2x2_2x2_2) {
  // example from Moler & Van Loan, 2003
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(2, 2), m2(2, 2);

  m1 << -49, 24, -64, 31;
  m2 << -.735759, .551819, -1.471518, 1.103638;

  expect_matrix_eq(m2, stan::math::matrix_exp_2x2(m1));
}
