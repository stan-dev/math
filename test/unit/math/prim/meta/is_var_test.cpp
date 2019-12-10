#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_var) {
  using stan::is_var;
  EXPECT_FALSE(is_var<int>::value);
  EXPECT_FALSE(is_var<double>::value);
}
