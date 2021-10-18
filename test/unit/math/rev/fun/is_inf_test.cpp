#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(AgradRev, is_inf) {
  using stan::math::is_inf;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  stan::math::var a(infinity);
  EXPECT_TRUE(is_inf(a));

  stan::math::var b(3.0);
  EXPECT_FALSE(is_inf(b));

  stan::math::var c(nan);
  EXPECT_FALSE(is_inf(c));
}
