#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(AgradMix, primitiveValueRevNested) {
  using stan::math::fvar;
  using stan::math::primitive_value;
  using stan::math::var;

  fvar<var> a = 5.2;
  EXPECT_FLOAT_EQ(5.2, primitive_value(a));

  // make sure all work together
  EXPECT_FLOAT_EQ(5.3, primitive_value(5.3));
  EXPECT_EQ(3, primitive_value(3));
}

TEST(AgradMix, primitiveValueNanRevNested) {
  using stan::math::fvar;
  using stan::math::primitive_value;
  using stan::math::var;
  double nan = std::numeric_limits<double>::quiet_NaN();

  fvar<var> a = nan;
  EXPECT_TRUE(stan::math::is_nan(primitive_value(a)));
  EXPECT_TRUE(stan::math::is_nan(primitive_value(nan)));
}
