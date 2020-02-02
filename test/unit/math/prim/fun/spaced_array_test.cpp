#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

void expect_linspaced_array(int K, double low, double high,
                            const std::vector<double>& expected) {
  std::vector<double> found = stan::math::linspaced_array(K, low, high);
  expect_std_vector_eq(expected, found);
}

TEST(MathFunctions, linspaced_array) {
  expect_linspaced_array(0, 1, 5, {});
  expect_linspaced_array(1, 1, 5, {5});
  expect_linspaced_array(5, 1, 5, {1, 2, 3, 4, 5});
  expect_linspaced_array(5, -2, 2, {-2, -1, 0, 1, 2});
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
