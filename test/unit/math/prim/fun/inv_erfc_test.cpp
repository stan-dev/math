#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, inv_erfc) {
  using stan::math::inv_erfc;
  EXPECT_FLOAT_EQ(-0.732869077959, inv_erfc(1.7));
  EXPECT_FLOAT_EQ(0, inv_erfc(1));
  EXPECT_FLOAT_EQ(1.82138636772, inv_erfc(0.01));
  EXPECT_FLOAT_EQ(-1.16308715368, inv_erfc(1.9));
}

TEST(MathFunctions, inv_erfcNan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::inv_erfc(nan)));
}
