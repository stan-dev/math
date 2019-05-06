#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, beta) {
  using stan::math::beta;

  EXPECT_FLOAT_EQ(beta(2.15, 1.71), 0.1936023967178879658641281697269);
  EXPECT_FLOAT_EQ(beta(7.62, 10.15), 0.0000065340071071564116445887286);
}

TEST(MathFunctions, beta_nan) {
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;
  using stan::math::beta;

  EXPECT_PRED1(boost::math::isnan<double>, beta(NOT_A_NUMBER, 2.16));
  EXPECT_PRED1(boost::math::isnan<double>, beta(1.65, NOT_A_NUMBER));

  EXPECT_PRED1(boost::math::isnan<double>, beta(INFTY, 2.16));
  EXPECT_PRED1(boost::math::isnan<double>, beta(1.65, INFTY));
}
