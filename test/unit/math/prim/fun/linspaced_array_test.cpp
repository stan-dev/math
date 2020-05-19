#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MathFunctions, linspaced_array) {
  using stan::math::linspaced_array;
  EXPECT_STD_VECTOR_FLOAT_EQ(linspaced_array(0, 1, 5), std::vector<int>({}));
  EXPECT_STD_VECTOR_FLOAT_EQ(linspaced_array(1, 1, 5), std::vector<int>({5}));
  EXPECT_STD_VECTOR_FLOAT_EQ(linspaced_array(5, 1, 5),
                             std::vector<int>({1, 2, 3, 4, 5}));
  EXPECT_STD_VECTOR_FLOAT_EQ(linspaced_array(5, -2, 2),
                             std::vector<int>({-2, -1, 0, 1, 2}));
}

TEST(MathFunctions, linspaced_array_throw) {
  using stan::math::linspaced_array;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  int K = 5;
  double low = -2;
  double high = 6;

  EXPECT_THROW(linspaced_array(-1, low, high), std::domain_error);

  EXPECT_THROW(linspaced_array(K, inf, high), std::domain_error);
  EXPECT_THROW(linspaced_array(K, nan, high), std::domain_error);

  EXPECT_THROW(linspaced_array(K, low, low - 1), std::domain_error);
  EXPECT_THROW(linspaced_array(K, low, inf), std::domain_error);
  EXPECT_THROW(linspaced_array(K, low, nan), std::domain_error);
}
