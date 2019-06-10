#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, length_mvt) {
  using stan::max_size;

  double x1 = 1.0, x2 = 0.0, x3 = 2.0, x4 = -3.0, x5 = -11.0;
  EXPECT_EQ(1U, max_size(x1, x2));
  EXPECT_EQ(1U, max_size(x1, x2, x3));
  EXPECT_EQ(1U, max_size(x1, x2, x3, x4));
  EXPECT_EQ(1U, max_size(x1, x2, x3, x4, x5));
}
