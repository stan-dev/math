#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>

TEST(MathFunctions, log1m) {
  EXPECT_FLOAT_EQ(stan::math::log1p(-0.1), stan::math::log1m(0.1));
}

TEST(MathFunctions, log1mOverflow) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(), stan::math::log1m(1.0));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(), stan::math::log1m(1));
}

TEST(MathFunctions, log1m_exception) {
  EXPECT_THROW_MSG(stan::math::log1m(10.0), std::domain_error,
                   "log1m: x is 10, but must be less than or equal to 1");
}

TEST(MathFunctions, log1m_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_PRED1(boost::math::isnan<double>, stan::math::log1m(nan));
}
