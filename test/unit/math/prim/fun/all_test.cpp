#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, all) {
  using stan::math::all;

  EXPECT_TRUE(all(true));
  EXPECT_TRUE(all(0 < 1));

  Eigen::ArrayXd inp(2);
  inp << 1, 2;

  EXPECT_FALSE(all(inp < 2));

  Eigen::Array<bool, -1, 1> bool_input(4);
  bool_input << true, true, true, true;

  EXPECT_TRUE(all(bool_input));

  bool_input(2) = false;

  EXPECT_FALSE(all(bool_input));

  bool_input(2) = true;

  std::vector<Eigen::Array<bool, -1, 1>> nested_bool(2);
  nested_bool[0] = bool_input;
  nested_bool[1] = bool_input;

  EXPECT_TRUE(all(nested_bool));

  nested_bool[1](3) = false;

  EXPECT_FALSE(all(nested_bool));
}
