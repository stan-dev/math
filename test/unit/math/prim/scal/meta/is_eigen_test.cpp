#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, primitive) {
  using stan::is_eigen;
  EXPECT_FALSE((is_eigen<bool>::value));
  EXPECT_FALSE((is_eigen<double>::value));
  EXPECT_FALSE((is_eigen<int>::value));
}
