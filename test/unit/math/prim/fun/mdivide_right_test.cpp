#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, mdivide_right_val) {
  stan::math::matrix_d Ad(2, 2);
  Ad << 2.0, 3.0, 5.0, 7.0;

  stan::math::matrix_d I = Eigen::MatrixXd::Identity(2, 2);
  EXPECT_MATRIX_FLOAT_EQ(I, stan::math::mdivide_left(Ad, Ad));
}

TEST(MathMatrixPrim, mdivide_right_val2) {
  stan::math::row_vector_d b(5);
  stan::math::matrix_d A(5, 5);
  stan::math::row_vector_d expected(5);

  b << 19, 150, -170, 140, 31;
  A << 1, 8, -9, 7, 5, 0, 1, 0, 4, 4, 0, 0, 1, 2, 5, 0, 0, 0, 1, -5, 0, 0, 0, 0,
      1;
  expected << 19, -2, 1, 13, 4;

  EXPECT_MATRIX_FLOAT_EQ(expected, stan::math::mdivide_right(b, A));
}

TEST(MathMatrixPrim, mdivide_right_size_zero) {
  using stan::math::mdivide_right;
  stan::math::matrix_d m1, m2, res;

  m1.resize(0, 2);
  m2.resize(2, 2);
  m2 << 3, 5, 7, 11;
  res = mdivide_right(m1, m2);
  EXPECT_EQ(m1.rows(), res.rows());
  EXPECT_EQ(m2.cols(), res.cols());

  m1.resize(2, 0);
  m2.resize(0, 0);
  res = mdivide_right(m1, m2);
  EXPECT_EQ(m1.rows(), res.rows());
  EXPECT_EQ(m2.cols(), res.cols());
}
