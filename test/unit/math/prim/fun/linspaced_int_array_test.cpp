#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MathFunctions, linspaced_int_array) {
  using stan::math::linspaced_int_array;
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(0, 1, 5), std::vector<int>({}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(1, 1, 5), std::vector<int>({5}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(5, 1, 5),
                       std::vector<int>({1, 2, 3, 4, 5}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(5, -2, 2),
                       std::vector<int>({-2, -1, 0, 1, 2}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(5, 0, 8),
                       std::vector<int>({0, 2, 4, 6, 8}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(5, 0, 8),
                       std::vector<int>({0, 2, 4, 6, 8}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(5, 0, 11),
                       std::vector<int>({0, 2, 4, 6, 8}));
  EXPECT_STD_VECTOR_EQ(linspaced_int_array(5, 0, 12),
                       std::vector<int>({0, 3, 6, 9, 12}));
}

TEST(MathFunctions, linspaced_int_array_throw) {
  using stan::math::linspaced_int_array;
  int low = -2;
  int high = 5;
  EXPECT_THROW(linspaced_int_array(-1, low, high), std::domain_error);
  EXPECT_THROW(linspaced_int_array(5, low, low - 1), std::domain_error);
}
