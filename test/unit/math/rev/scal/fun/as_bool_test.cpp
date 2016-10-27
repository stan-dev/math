#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>

TEST(AgradRev,asBool) {
  using stan::math::as_bool;
  using stan::math::var;

  EXPECT_TRUE(as_bool(var(1)));
  EXPECT_TRUE(as_bool(var(-10L)));
  EXPECT_TRUE(as_bool(var(1.7)));
  EXPECT_TRUE(as_bool(var(-1.7)));
  EXPECT_TRUE(as_bool(var(std::numeric_limits<double>::infinity())));
  EXPECT_TRUE(as_bool(var(-std::numeric_limits<double>::infinity())));
  // don't like this behavior, but it's what C++ does
  EXPECT_TRUE(as_bool(var(std::numeric_limits<double>::quiet_NaN())));

  EXPECT_FALSE(as_bool(var(0)));
  EXPECT_FALSE(as_bool(var(0.0)));
  EXPECT_FALSE(as_bool(var(0.0f)));
}
TEST(AgradRev,as_bool_nan) {
  stan::math::var nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(stan::math::as_bool(nan));
}
