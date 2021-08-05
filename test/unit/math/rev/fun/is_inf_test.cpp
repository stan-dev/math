#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(AgradRev, is_inf) {
  using stan::math::is_inf;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  AVAR a(infinity);
  EXPECT_TRUE(is_inf(a));

  AVAR b(3.0);
  EXPECT_FALSE(is_inf(b));

  AVAR c(nan);
  EXPECT_FALSE(is_inf(c));
}
