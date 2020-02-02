#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <limits>

void expect_range_row_vector(double low, double high, double by,
                             const Eigen::RowVectorXd& expected) {
  Eigen::RowVectorXd found = stan::math::range_row_vector(low, high, by);
  expect_matrix_eq(expected, found);
}

TEST(MathFunctions, range_row_vector) {
  Eigen::RowVectorXd v(1);
  v << 1;
  expect_range_row_vector(1, 1, 1, v);
  expect_range_row_vector(1, 2, 1, v);
  expect_range_row_vector(1, 2, 5, v);

  double by = 0.5;
  Eigen::RowVectorXd v1(5);
  v1 << 1, 1.5, 2, 2.5, 3;
  expect_range_row_vector(1, 3.49, by, v1);
  expect_range_row_vector(1, 3.50, by, v1);

  Eigen::RowVectorXd v2(6);
  v2 << v1, 3.5;
  expect_range_row_vector(1, 3.51, by, v2);
}

TEST(MathFunctions, range_row_vector_throw) {
  using stan::math::range_row_vector;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double by = 0.5;
  double low = -2;
  double high = 6;

  EXPECT_THROW(range_row_vector(low, high, 0), std::domain_error);
  EXPECT_THROW(range_row_vector(low, high, -1), std::domain_error);

  EXPECT_THROW(range_row_vector(inf, high, by), std::domain_error);
  EXPECT_THROW(range_row_vector(nan, high, by), std::domain_error);

  EXPECT_THROW(range_row_vector(low, low - 1, by), std::domain_error);
  EXPECT_THROW(range_row_vector(low, inf, by), std::domain_error);
  EXPECT_THROW(range_row_vector(low, nan, by), std::domain_error);
}
