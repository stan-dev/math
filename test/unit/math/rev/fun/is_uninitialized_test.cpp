#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST_F(AgradRev, Rev_undefined) {
  stan::math::var a;
  EXPECT_TRUE(a.is_uninitialized());
  a = 5;
  EXPECT_FALSE(a.is_uninitialized());
}

TEST_F(AgradRev, Rev_is_uninitialized_nan) {
  stan::math::var nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FALSE(stan::math::is_uninitialized(nan));
}
