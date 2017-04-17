#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/util.hpp>

TEST(MathFunctions, sign) {
  using stan::math::var;
  var x;
  x = 0;
  EXPECT_EQ(0, stan::math::sign(x));
  x = 0.0000001;
  EXPECT_EQ(1, stan::math::sign(x));
  x = -0.001;
  EXPECT_EQ(-1, stan::math::sign(x));
}

TEST(MathFunctions, check_varis_on_stack) {
  stan::math::var x = 0;
  test::check_varis_on_stack(stan::math::sign(x));
}
