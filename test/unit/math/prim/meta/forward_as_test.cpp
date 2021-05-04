#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, ForwardAsScalar) {
  using stan::math::forward_as;

  int a = 1;

  EXPECT_NO_THROW(forward_as<int>(a) = 2);
  EXPECT_EQ(a, 2);
  EXPECT_THROW(forward_as<double>(a), std::runtime_error);
}
