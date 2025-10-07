#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST_F(AgradRev, Rev_is_nan) {
  using stan::math::is_nan;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  stan::math::var a(nan);
  EXPECT_TRUE(is_nan(a));

  stan::math::var b(3.0);
  EXPECT_FALSE(is_nan(b));

  stan::math::var c(infinity);
  EXPECT_FALSE(is_nan(c));
}
