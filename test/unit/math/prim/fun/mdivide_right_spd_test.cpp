#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, mdivide_right_spd_val) {
  using stan::math::mdivide_right_spd;
  stan::math::matrix_d Ad(2, 2);
  Ad << 2.0, 3.0, 3.0, 7.0;

  stan::math::matrix_d I = mdivide_right_spd(Ad, Ad);
  EXPECT_NEAR(1.0, I(0, 0), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1), 1.0e-12);

  Ad << 1, 1, 1, 1;
  EXPECT_THROW(mdivide_right_spd(Ad, Ad), std::domain_error);
}

TEST(MathMatrixPrim, mdivide_right_spd_size_zero) {
  using stan::math::mdivide_right_spd;
  stan::math::matrix_d m1, m2, res;

  m1.resize(0, 2);
  m2.resize(2, 2);
  m2 << 7, 2, 2, 4;
  res = mdivide_right_spd(m1, m2);
  EXPECT_EQ(m1.rows(), res.rows());
  EXPECT_EQ(m2.cols(), res.cols());

  m1.resize(2, 0);
  m2.resize(0, 0);
  res = mdivide_right_spd(m1, m2);
  EXPECT_EQ(m1.rows(), res.rows());
  EXPECT_EQ(m2.cols(), res.cols());
}
