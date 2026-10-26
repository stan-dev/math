#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, is_even_scalar) {
  EXPECT_EQ(1, stan::math::is_even(0));
  EXPECT_EQ(0, stan::math::is_even(1));
  EXPECT_EQ(1, stan::math::is_even(2));
  EXPECT_EQ(0, stan::math::is_even(3));

  EXPECT_EQ(0, stan::math::is_even(-1));
  EXPECT_EQ(1, stan::math::is_even(-2));
  EXPECT_EQ(0, stan::math::is_even(-3));
}

TEST(MathFunctions, is_even_vector) {
  std::vector<int> input{0, 1, 2, 3, -3, -4};
  std::vector<int> expected{1, 0, 1, 0, 0, 1};

  EXPECT_EQ(expected, stan::math::is_even(input));
}

TEST(MathFunctions, is_even_nested_vector) {
  std::vector<std::vector<int>> input{{1, 2}, {3, 4}};

  std::vector<std::vector<int>> expected{{0, 1}, {0, 1}};

  EXPECT_EQ(expected, stan::math::is_even(input));
}

TEST(MathFunctions, is_even_empty_vector) {
  std::vector<int> input;
  std::vector<int> expected;

  EXPECT_EQ(expected, stan::math::is_even(input));
}
