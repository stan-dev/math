#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, length) {
  using stan::length;
  std::vector<double> x(10);
  EXPECT_EQ(10U, length(x));
}
