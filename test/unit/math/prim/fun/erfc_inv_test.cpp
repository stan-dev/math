#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, erfc_inv) {
  using stan::math::erfc_inv;
  EXPECT_FLOAT_EQ(-0.732869077959, erfc_inv(1.7));
  EXPECT_FLOAT_EQ(0, erfc_inv(1));
  EXPECT_FLOAT_EQ(1.82138636772, erfc_inv(0.01));
  EXPECT_FLOAT_EQ(-1.16308715368, erfc_inv(1.9));
}

TEST(MathFunctions, erfc_invNan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::erfc_inv(nan)));
}
