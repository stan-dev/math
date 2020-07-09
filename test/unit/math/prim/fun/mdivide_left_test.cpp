#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, mdivide_left_val) {
  stan::math::matrix_d Ad(2, 2);
  Ad << 2.0, 3.0, 5.0, 7.0;

  stan::math::matrix_d I = Eigen::MatrixXd::Identity(2, 2);
  EXPECT_MATRIX_FLOAT_EQ(I, stan::math::mdivide_left(Ad, Ad));
}

TEST(MathMatrixPrim, mdivide_left_val2) {
  stan::math::matrix_d A(5, 5);
  stan::math::vector_d b(5);
  stan::math::vector_d expected(5);

  A << 1, 8, -9, 7, 5, 0, 1, 0, 4, 4, 0, 0, 1, 2, 5, 0, 0, 0, 1, -5, 0, 0, 0, 0,
      1;
  b << 19, 150, -170, 140, 31;
  expected << -1204, -1154, -915, 295, 31;

  EXPECT_MATRIX_FLOAT_EQ(expected, stan::math::mdivide_left(A, b));
}

TEST(MathMatrixPrim, mdivide_left_size_zero) {
  using stan::math::mdivide_left;
  stan::math::matrix_d m1, m2, res;

  m1.resize(2, 2);
  m1 << 3, 5, 7, 11;
  m2.resize(2, 0);
  res = mdivide_left(m1, m2);
  EXPECT_EQ(m1.rows(), res.rows());
  EXPECT_EQ(m2.cols(), res.cols());

  m1.resize(0, 0);
  m2.resize(0, 2);
  res = mdivide_left(m1, m2);
  EXPECT_EQ(m1.rows(), res.rows());
  EXPECT_EQ(m2.cols(), res.cols());
}
