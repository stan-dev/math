#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, acosh) {
  using stan::math::acosh;
  EXPECT_FLOAT_EQ(0, acosh(1));
  EXPECT_FLOAT_EQ(0.96242365, acosh(1.5));
  EXPECT_FLOAT_EQ(3.0797991, acosh(10.9));
}

TEST(MathFunctions, acosh_exception) {
  using stan::math::acosh;
  EXPECT_THROW(acosh(0.5), std::domain_error);
}

TEST(MathFunctions, acosh_inf_return) {
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::acosh(std::numeric_limits<double>::infinity()));
}

TEST(MathFunctions, acosh_nan) {
  using stan::math::acosh;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::acosh(nan)));
}
