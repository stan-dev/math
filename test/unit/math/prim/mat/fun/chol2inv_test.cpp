#include <stan/math/prim/mat.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, chol2inv_exception) {
  using stan::math::chol2inv;

  stan::math::matrix_d m1(2, 3);

  // non-square
  m1 << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(chol2inv(m1), std::invalid_argument);

  stan::math::matrix_d m2(3, 3);

  // non-lower-triangular
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  EXPECT_THROW(chol2inv(m2), std::domain_error);
}

TEST(MathMatrix, chol2inv) {
  using stan::math::chol2inv;
  using stan::math::cholesky_decompose;
  using stan::math::inverse_spd;
  using stan::math::matrix_d;
  using stan::math::wishart_rng;

  boost::random::mt19937 rng;
  matrix_d I(3, 3);
  I.setZero();
  I.diagonal().setOnes();
  matrix_d Y = wishart_rng(4.0, I, rng);
  matrix_d L = cholesky_decompose(Y);
  matrix_d Y_inv = inverse_spd(Y);
  matrix_d Y_inv2 = chol2inv(L);
  for (int j = 0; j < Y.cols(); j++)
    for (int i = 0; i < Y.rows(); i++)
      EXPECT_FLOAT_EQ(Y_inv(i, j), Y_inv2(i, j));
}

TEST(MathMatrix, chol2inv01) {
  using stan::math::chol2inv;
  using stan::math::matrix_d;

  matrix_d Y(0, 0);
  matrix_d Y_inv2 = chol2inv(Y);
  EXPECT_EQ(0, Y_inv2.rows());
  EXPECT_EQ(0, Y_inv2.cols());

  matrix_d L(1, 1);
  L(0, 0) = 3.0;
  matrix_d inv2 = chol2inv(L);
  EXPECT_FLOAT_EQ(1 / 9.0, inv2(0, 0));
}
