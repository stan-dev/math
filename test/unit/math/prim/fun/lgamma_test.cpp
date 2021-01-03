#include <stan/math/prim.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, lgamma) { EXPECT_TRUE(stan::math::is_inf(lgamma(0.0))); }

TEST(MathFunctions, lgammaStanMathUsing) { using stan::math::lgamma; }

TEST(MathFunctions, lgamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::lgamma(nan)));
  EXPECT_TRUE(std::isinf(stan::math::lgamma(0)));
  // up to boost 1.70.0 the boost::math::lgamma return NaN for large
  // arguments (large is greater than sqrt(largest double of 1E308)
  // when used with the stan::math::boost_policy_t which avoids
  // propagation of input arguments from double to long double
  // internally to boost::math::lgamma. See
  // https://github.com/boostorg/math/issues/242 for this. The
  // stan::math::lgamma implementation is based on the
  // boost::math::lgamma only when MinGW compilers are used. Thus, to
  // ensure that boost::math::lgamma contains the needed bugfixes we
  // test here specifically the boost::math::lgamma by testing for a
  // finite return for a large argument.
  EXPECT_TRUE(std::isnormal(
      boost::math::lgamma(1.0E50, stan::math::boost_policy_t<>())));
}
