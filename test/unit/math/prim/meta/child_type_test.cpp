#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, value_type) {
  using stan::math::child_type;

  EXPECT_SAME_TYPE(double, child_type<double>::type);
}
