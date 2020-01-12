#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::diag_post_multiply;

TEST(MathMatrixPrimMat, diagPostMultiply) {
  Matrix<double, Dynamic, Dynamic> m(1, 1);
  m << 3;

  Matrix<double, Dynamic, 1> v(1);
  v << 9;

  Matrix<double, Dynamic, Dynamic> v_m(1, 1);
  v_m << 9;

  expect_matrix_eq(m * v_m, diag_post_multiply(m, v));
}
TEST(MathMatrixPrimMat, diagPostMultiply2) {
  Matrix<double, Dynamic, Dynamic> m(2, 2);
  m << 2, 3, 4, 5;

  Matrix<double, Dynamic, 1> v(2);
  v << 10, 100;

  Matrix<double, Dynamic, Dynamic> v_m(2, 2);
  v_m << 10, 0, 0, 100;

  expect_matrix_eq(m * v_m, diag_post_multiply(m, v));

  Matrix<double, 1, Dynamic> rv(2);
  rv << 10, 100;
  expect_matrix_eq(m * v_m, diag_post_multiply(m, rv));
}

TEST(MathMatrixPrimMat, diagPostMultiplyException) {
  Matrix<double, Dynamic, Dynamic> m(2, 2);
  m << 2, 3, 4, 5;
  EXPECT_THROW(diag_post_multiply(m, m), std::invalid_argument);

  Matrix<double, Dynamic, 1> v(3);
  v << 1, 2, 3;
  EXPECT_THROW(diag_post_multiply(m, v), std::invalid_argument);
}
