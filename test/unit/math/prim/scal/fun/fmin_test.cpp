#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, fminFinite) {
  using stan::math::fmin;
  EXPECT_FLOAT_EQ(1.0, fmin(1, 2));
  EXPECT_FLOAT_EQ(1.0, fmin(1.0, 2));
  EXPECT_FLOAT_EQ(1.0, fmin(1, 2.0));
  EXPECT_FLOAT_EQ(1.0, fmin(1.0, 2.0));

  EXPECT_FLOAT_EQ(1.0, fmin(2, 1));
  EXPECT_FLOAT_EQ(1.0, fmin(2, 1.0));
  EXPECT_FLOAT_EQ(1.0, fmin(2.0, 1));
  EXPECT_FLOAT_EQ(1.0, fmin(2.0, 1.0));
}

TEST(MathFunctions, fminNaN) {
  using stan::math::fmin;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FLOAT_EQ(1.0, fmin(1, nan));
  EXPECT_FLOAT_EQ(1.0, fmin(nan, 1));
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::fmin(nan, nan));
}

TEST(MathFunctions, fminInf) {
  using stan::math::fmin;
  double inf = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1, fmin(inf, 1));
  EXPECT_FLOAT_EQ(1, fmin(1, inf));
  EXPECT_FLOAT_EQ(-inf, fmin(inf, -inf));
  EXPECT_FLOAT_EQ(-inf, fmin(-inf, inf));
}
