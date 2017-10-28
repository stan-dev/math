#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

void test_logical_or(double x, double y) {
  AVAR x_v = x;
  AVAR y_v = y;
  EXPECT_EQ(x || y, x_v || y_v);
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
      test_logical_or(xs[i], xs[j]);
}
