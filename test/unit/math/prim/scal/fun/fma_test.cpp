#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

// this is just testing the nan behavior of the built-in fma
// there is no longer a stan::math::fma, just the agrad versions
// instead, the top-level ::fma should be used by including <cmath>

TEST(MathFunctions, fma) {
  using stan::math::fma;
  EXPECT_FLOAT_EQ(5.0, fma(1.0, 2.0, 3.0));
  EXPECT_FLOAT_EQ(10.0, fma(2.0, 3.0, 4.0));
  EXPECT_FLOAT_EQ(
      11.0, fma(static_cast<int>(3), static_cast<int>(2), static_cast<int>(5)));
}

TEST(MathFunctions, fma_nan) {
  using stan::math::fma;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(fma(1.0, 2.0, nan)));

  EXPECT_TRUE(std::isnan(fma(1.0, nan, 3.0)));

  EXPECT_TRUE(std::isnan(fma(1.0, nan, nan)));

  EXPECT_TRUE(std::isnan(fma(nan, 2.0, 3.0)));

  EXPECT_TRUE(std::isnan(fma(nan, 2.0, nan)));

  EXPECT_TRUE(std::isnan(fma(nan, nan, 3.0)));

  EXPECT_TRUE(std::isnan(fma(nan, nan, nan)));
}
