#include <stan/math/rev/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaRevScal, value_type) {
  using stan::math::child_type;
  using stan::math::var;

  EXPECT_SAME_TYPE(double, child_type<var>::type);
}
