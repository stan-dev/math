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

  Eigen::Array<bool, -1, 1> bool_array(4);
  bool_array << false, false, false, false;

  std::vector<bool> bool_stdvector{false, false};

  std::vector<Eigen::Array<bool, -1, 1>> nested_bool(2);
  nested_bool[0] = bool_array;
  nested_bool[1] = bool_array;

  auto bool_tuple
      = std::make_tuple(false, bool_stdvector, bool_array, nested_bool);
  std::vector<decltype(bool_tuple)> stdvec_bool_tuple{bool_tuple};

  EXPECT_TRUE(any(inp < 2));
  EXPECT_FALSE(any(bool_stdvector));
  EXPECT_FALSE(any(bool_array));
  EXPECT_FALSE(any(nested_bool));
  EXPECT_FALSE(any(bool_tuple));
  EXPECT_FALSE(any(stdvec_bool_tuple));

  bool_array(2) = true;
  nested_bool[1](3) = true;
  bool_stdvector[1] = true;
  std::get<3>(bool_tuple) = nested_bool;
  stdvec_bool_tuple[0] = bool_tuple;

  EXPECT_TRUE(any(bool_array));
  EXPECT_TRUE(any(nested_bool));
  EXPECT_TRUE(any(bool_stdvector));
  EXPECT_TRUE(any(bool_tuple));
  EXPECT_TRUE(any(stdvec_bool_tuple));
}
