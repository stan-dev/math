#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, size_of) {
  using stan::size_of;
  double x1 = 2;
  EXPECT_EQ(1U, size_of(x1));
}
