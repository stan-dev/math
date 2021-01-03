#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

void expect_linspaced_vector(int K, double low, double high,
                             const Eigen::VectorXd& expected) {
  Eigen::VectorXd found = stan::math::linspaced_vector(K, low, high);
  EXPECT_MATRIX_FLOAT_EQ(expected, found);
}

TEST(MathFunctions, linspaced_vector) {
  expect_linspaced_vector(0, 1, 5, {});

  Eigen::VectorXd v(1);
  v << 5;
  expect_linspaced_vector(1, 1, 5, v);

  int K = 5;
  Eigen::VectorXd v1(K);
  v1 << 1, 2, 3, 4, 5;
  expect_linspaced_vector(K, 1, 5, v1);

  Eigen::VectorXd v2(K);
  v2 << -2, -1, 0, 1, 2;
  expect_linspaced_vector(K, -2, 2, v2);
}

TEST(MathFunctions, linspaced_vector_throw) {
  using stan::math::linspaced_vector;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  int K = 5;
  double low = -2;
  double high = 6;

  EXPECT_THROW(linspaced_vector(-1, low, high), std::domain_error);

  EXPECT_THROW(linspaced_vector(K, inf, high), std::domain_error);
  EXPECT_THROW(linspaced_vector(K, nan, high), std::domain_error);

  EXPECT_THROW(linspaced_vector(K, low, low - 1), std::domain_error);
  EXPECT_THROW(linspaced_vector(K, low, inf), std::domain_error);
  EXPECT_THROW(linspaced_vector(K, low, nan), std::domain_error);
}
