#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, set_spaced_row_vector) {
  using Eigen::RowVectorXd;
  using stan::math::set_spaced_row_vector;

  RowVectorXd u0 = set_spaced_row_vector(0, 1, 2);
  EXPECT_EQ(0, u0.size());

  double low = 1;
  double high = 5;
  RowVectorXd u11 = set_spaced_row_vector(1, low, high);
  EXPECT_EQ(1, u11.size());
  EXPECT_FLOAT_EQ(high, u11[0]);

  int K = 5;
  RowVectorXd u1 = set_spaced_row_vector(K, 1, 5);
  RowVectorXd v1(K);
  v1 << 1, 2, 3, 4, 5;

  EXPECT_EQ(K, u1.size());
  for (int i = 0; i < K; i++) {
    EXPECT_FLOAT_EQ(u1[i], v1[i]);
  }

  RowVectorXd u2 = set_spaced_row_vector(K, -2, 2);
  RowVectorXd v2(K);
  v2 << -2, -1, 0, 1, 2;

  EXPECT_EQ(K, u2.size());
  for (int i = 0; i < K; i++) {
    EXPECT_FLOAT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, set_spaced_row_vector_throw) {
  using stan::math::set_spaced_row_vector;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  int K = 5;
  double low = -2;
  double high = 6;

  EXPECT_THROW(set_spaced_row_vector(-1, low, high), std::domain_error);

  EXPECT_THROW(set_spaced_row_vector(K, inf, high), std::domain_error);
  EXPECT_THROW(set_spaced_row_vector(K, nan, high), std::domain_error);

  EXPECT_THROW(set_spaced_row_vector(K, low, low - 1), std::domain_error);
  EXPECT_THROW(set_spaced_row_vector(K, low, inf), std::domain_error);
  EXPECT_THROW(set_spaced_row_vector(K, low, nan), std::domain_error);
}
