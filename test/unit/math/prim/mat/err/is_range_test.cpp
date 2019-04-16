#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(ErrorHandlingMatrix, isRange_6_arg_std_vector) {
  std::vector<double> x;

  x.resize(4);
  EXPECT_TRUE(stan::math::is_range(4, 1));
  EXPECT_TRUE(stan::math::is_range(4, 2));
  EXPECT_TRUE(stan::math::is_range(4, 3));
  EXPECT_TRUE(stan::math::is_range(4, 4));
  EXPECT_FALSE(stan::math::is_range(4, 12));
}

TEST(ErrorHandlingMatrix, isRange_4_arg_std_vector) {
  std::vector<double> x;

  x.resize(4);
  EXPECT_TRUE(stan::math::is_range(4, 1));
  EXPECT_TRUE(stan::math::is_range(4, 2));
  EXPECT_TRUE(stan::math::is_range(4, 3));
  EXPECT_TRUE(stan::math::is_range(4, 4));
  EXPECT_FALSE(stan::math::is_range(4, 12));
}
