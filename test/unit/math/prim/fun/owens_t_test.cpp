#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, owens_t) {
  double a = 1.0;
  double b = 2.0;
  EXPECT_FLOAT_EQ(stan::math::owens_t(a, b), boost::math::owens_t(a, b));
}

TEST(MathFunctions, owens_t_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::owens_t(1.0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::owens_t(nan, 2.0)));
  EXPECT_TRUE(std::isnan(stan::math::owens_t(nan, nan)));
}
