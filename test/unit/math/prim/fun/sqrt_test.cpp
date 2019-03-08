#include <stan/math/prim.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, sqrtInt) {
  using stan::math::sqrt;
  using std::sqrt;
  EXPECT_FLOAT_EQ(std::sqrt(3.0), sqrt(3));
  EXPECT_TRUE(stan::math::is_nan(sqrt(-2)));
}
