#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, any) {
  using stan::math::any;

  EXPECT_TRUE(any(true));
  EXPECT_TRUE(any(0 < 1));

  Eigen::ArrayXd inp(2);
  inp << 1, 2;

  EXPECT_TRUE(any(inp < 2));

  Eigen::Array<bool, -1, 1> bool_input(4);
  bool_input << false, false, false, false;

  EXPECT_FALSE(any(bool_input));

  bool_input(2) = true;

  EXPECT_TRUE(any(bool_input));

  bool_input(2) = false;

  std::vector<Eigen::Array<bool, -1, 1>> nested_bool(2);
  nested_bool[0] = bool_input;
  nested_bool[1] = bool_input;

  EXPECT_FALSE(any(nested_bool));

  nested_bool[1](3) = true;

  EXPECT_TRUE(any(nested_bool));
}
