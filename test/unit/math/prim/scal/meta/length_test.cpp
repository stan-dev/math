#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::length;

TEST(MetaTraitsPrimScal, length) {
  using stan::length;
  EXPECT_EQ(1U, length(27.0));
  EXPECT_EQ(1U, length(3));
}

TEST(MetaTraitsPrimArr, length) {
  using stan::length;
  std::vector<double> x(10);
  EXPECT_EQ(10U, length(x));
}
