#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, lgamma) { EXPECT_TRUE(stan::math::is_inf(lgamma(0.0))); }

TEST(MathFunctions, lgammaStanMathUsing) { using stan::math::lgamma; }

TEST(MathFunctions, lgamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::lgamma(nan));
  EXPECT_PRED1(boost::math::isinf<double>, stan::math::lgamma(0));
}
