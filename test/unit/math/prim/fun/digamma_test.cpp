#include <stan/math/prim.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, digamma) {
  EXPECT_FLOAT_EQ(boost::math::digamma(0.5), stan::math::digamma(0.5));
  EXPECT_FLOAT_EQ(boost::math::digamma(-1.5), stan::math::digamma(-1.5));
}

TEST(MathFunctions, digamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::digamma(nan)));

  EXPECT_TRUE(std::isnan(stan::math::digamma(-1)));

  EXPECT_TRUE(std::isnormal(stan::math::digamma(1.0E50)));
}

TEST(MathFunctions, digamma_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::digamma(b));
}
