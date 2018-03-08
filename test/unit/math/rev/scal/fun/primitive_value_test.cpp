#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

TEST(AgradRev, primitiveValue) {
  using stan::math::primitive_value;
  using stan::math::var;

  var a = 5.0;
  EXPECT_FLOAT_EQ(5.0, primitive_value(a));

  // make sure all work together
  EXPECT_FLOAT_EQ(5.0, primitive_value(5.0));
  EXPECT_EQ(5, primitive_value(5));
}
