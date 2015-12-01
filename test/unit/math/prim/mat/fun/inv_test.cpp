#include <gtest/gtest.h>
#include <stan/math/prim/mat/fun/inv.hpp>
#include <vector>

TEST(primMatFun, inv) {
  using stan::math::inv;
  std::vector<double> x;
  x.push_back(2);
  x.push_back(4);
  std::vector<double> y = inv(x);
  EXPECT_EQ(x.size(), y.size());
  EXPECT_FLOAT_EQ(0.5, y[0]);
  EXPECT_FLOAT_EQ(0.25, y[1]);
  
}
