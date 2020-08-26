#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST_F(AgradRev, is_nan) {
  using stan::math::is_nan;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  AVAR a(nan);
  EXPECT_TRUE(is_nan(a));

  AVAR b(3.0);
  EXPECT_FALSE(is_nan(b));

  AVAR c(infinity);
  EXPECT_FALSE(is_nan(c));
}
