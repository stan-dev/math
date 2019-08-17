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
  // up to boost 1.70.0 the boost::math::lgamma return NaN for large
  // arguments (large is creater than sqrt(largest double) when used
  // with the policy which avoids propagation of double to long double
  // see https://github.com/boostorg/math/issues/242
  // this boost::math::lgamma is used and as such we need to make sure
  // that the appropiate boost is available
  EXPECT_PRED1(boost::math::isnormal<double>,
               boost::math::lgamma(1.0E50, stan::math::boost_policy_t()));
}
