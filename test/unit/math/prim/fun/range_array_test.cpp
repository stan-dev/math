#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

void expect_range_array(double low, double high, double by,
                        const std::vector<double>& expected) {
  std::vector<double> found = stan::math::range_array(low, high, by);
  expect_std_vector_eq(expected, found);
}

TEST(MathFunctions, range_array) {
  std::vector<double> v{1};
  expect_range_array(1, 1, 1, v);
  expect_range_array(1, 2, 1, v);
  expect_range_array(1, 2, 5, v);

  double by = 0.5;
  std::vector<double> v1{1, 1.5, 2, 2.5, 3};
  expect_range_array(1, 3.49, by, v1);
  expect_range_array(1, 3.50, by, v1);

  std::vector<double> v2{1, 1.5, 2, 2.5, 3, 3.5};
  expect_range_array(1, 3.51, by, v2);
}

TEST(MathFunctions, range_array_throw) {
  using stan::math::range_array;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double by = 0.5;
  double low = -2;
  double high = 6;

  EXPECT_THROW(range_array(low, high, 0), std::domain_error);
  EXPECT_THROW(range_array(low, high, -1), std::domain_error);

  EXPECT_THROW(range_array(inf, high, by), std::domain_error);
  EXPECT_THROW(range_array(nan, high, by), std::domain_error);

  EXPECT_THROW(range_array(low, low - 1, by), std::domain_error);
  EXPECT_THROW(range_array(low, inf, by), std::domain_error);
  EXPECT_THROW(range_array(low, nan, by), std::domain_error);
}
