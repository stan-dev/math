#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, is_fvar) {
  using stan::is_fvar;
  EXPECT_FALSE(is_fvar<int>::value);
  EXPECT_FALSE(is_fvar<double>::value);
}
