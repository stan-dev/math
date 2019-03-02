#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

using stan::math::is_ordered;

TEST(ErrorHandling, isOrdered) {
  std::vector<double> y_;
  y_.push_back(0.0);
  y_.push_back(1.0);
  y_.push_back(2.0);
  EXPECT_TRUE(is_ordered(y_));

  y_[2] = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_ordered(y_));

  y_[0] = -10.0;
  EXPECT_TRUE(is_ordered(y_));

  y_[0] = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(is_ordered(y_));

  y_[0] = 0.0;
  y_[1] = 0.0;
  y_[2] = 0.0;
  EXPECT_FALSE(is_ordered(y_));

  y_[1] = std::numeric_limits<double>::infinity();
  y_[2] = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(is_ordered(y_));
}

TEST(ErrorHandling, isOrdered_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> y_;
  y_.push_back(0.0);
  y_.push_back(1.0);
  y_.push_back(2.0);
  for (size_t i = 0; i < y_.size(); i++) {
    y_[i] = nan;
    EXPECT_FALSE(is_ordered(y_));
    y_[i] = i;
  }
  for (size_t i = 0; i < y_.size(); i++) {
    y_[0] = 0.0;
    y_[1] = 10.0;
    y_[2] = std::numeric_limits<double>::infinity();
    y_[i] = nan;
    EXPECT_FALSE(is_ordered(y_));
  }
}
