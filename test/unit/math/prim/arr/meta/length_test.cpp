#include <gtest/gtest.h>
#include <stan/math/prim/arr.hpp>
#include <vector>

TEST(MetaTraits, length) {
  using stan::length;
  std::vector<double> x(10);
  EXPECT_EQ(10U, length(x));
}
