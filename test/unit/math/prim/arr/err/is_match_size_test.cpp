#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandling, isMatchingSize) {
  std::vector<double> a;
  std::vector<double> b;

  EXPECT_TRUE(stan::math::is_matching_size(a, b));

  a.push_back(3.0);
  a.push_back(3.0);

  EXPECT_FALSE(stan::math::is_matching_size(a, b));

  b.push_back(3.0);
  b.push_back(3.0);
  EXPECT_TRUE(stan::math::is_matching_size(a, b));
}

TEST(ErrorHandling, isMatchingSize_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> a;
  std::vector<double> b;
  EXPECT_TRUE(stan::math::is_matching_size(a, b));

  a.push_back(nan);
  a.push_back(nan);
  EXPECT_FALSE(stan::math::is_matching_size(a, b));

  b.push_back(nan);
  b.push_back(nan);
  EXPECT_TRUE(stan::math::is_matching_size(a, b));
}
