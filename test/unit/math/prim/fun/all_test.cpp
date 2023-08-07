#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, all) {
  using stan::math::all;

  EXPECT_TRUE(all(true));
  EXPECT_TRUE(all(0 < 1));

  Eigen::ArrayXd inp(2);
  inp << 1, 2;

  Eigen::Array<bool, -1, 1> bool_array(4);
  bool_array << true, true, true, true;

  std::vector<bool> bool_stdvector{true, true};

  std::vector<Eigen::Array<bool, -1, 1>> nested_bool(2);
  nested_bool[0] = bool_array;
  nested_bool[1] = bool_array;

  auto bool_tuple
      = std::make_tuple(true, bool_stdvector, bool_array, nested_bool);
  std::vector<decltype(bool_tuple)> stdvec_bool_tuple{bool_tuple};

  EXPECT_FALSE(all(inp < 2));
  EXPECT_TRUE(all(bool_stdvector));
  EXPECT_TRUE(all(bool_array));
  EXPECT_TRUE(all(nested_bool));
  EXPECT_TRUE(all(bool_tuple));
  EXPECT_TRUE(all(stdvec_bool_tuple));

  bool_array(2) = false;
  nested_bool[1](3) = false;
  bool_stdvector[1] = false;
  std::get<3>(bool_tuple) = nested_bool;
  stdvec_bool_tuple[0] = bool_tuple;

  EXPECT_FALSE(all(bool_array));
  EXPECT_FALSE(all(nested_bool));
  EXPECT_FALSE(all(bool_stdvector));
  EXPECT_FALSE(all(bool_tuple));
  EXPECT_FALSE(all(stdvec_bool_tuple));
}
