#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(ErrorHandling, isOrdered) {
  using stan::math::is_ordered;
  std::vector<double> y = {0, 1, 2};
  EXPECT_TRUE(is_ordered(y));

  y = {0, 1, std::numeric_limits<double>::infinity()};
  EXPECT_TRUE(is_ordered(y));

  y = {-10, 1, std::numeric_limits<double>::infinity()};
  EXPECT_TRUE(is_ordered(y));

  y = {-std::numeric_limits<double>::infinity(), 1,
       std::numeric_limits<double>::infinity()};
  EXPECT_TRUE(is_ordered(y));

  y = {0, 0, 0};
  EXPECT_FALSE(is_ordered(y));

  y = {0, std::numeric_limits<double>::infinity(),
       std::numeric_limits<double>::infinity()};
  EXPECT_FALSE(is_ordered(y));
}

TEST(ErrorHandling, isOrdered_nan) {
  using stan::math::is_ordered;
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> y = {0, 1, 2};

  for (size_t i = 0; i < y.size(); i++) {
    y[i] = nan;
    EXPECT_FALSE(is_ordered(y));
    y[i] = i;
  }
  for (size_t i = 0; i < y.size(); i++) {
    y = {0.0, 10.0, std::numeric_limits<double>::infinity()};
    y[i] = nan;
    EXPECT_FALSE(is_ordered(y));
  }
}
