#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

void test_logical_and(double x, double y) {
  using stan::math::var;
  EXPECT_EQ(x && y, var(x) && var(y));
  EXPECT_EQ(x && y, x && var(y));
  EXPECT_EQ(x && y, var(x) && y);
}

TEST(AgradRev, unaryNot) {
  std::vector<double> xs;
  xs.push_back(6.1);
  xs.push_back(6.1);
  xs.push_back(0);
  xs.push_back(-13.2);
  xs.push_back(std::numeric_limits<double>::infinity());
  xs.push_back(-std::numeric_limits<double>::infinity());
  xs.push_back(std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < xs.size(); ++i)
    for (size_t j = 0; j < xs.size(); ++j)
      test_logical_and(xs[i], xs[j]);
}
