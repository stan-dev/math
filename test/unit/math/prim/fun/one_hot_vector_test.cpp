#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, one_hot_vector) {
  using Eigen::VectorXd;
  using stan::math::one_hot_vector;

  int K = 4;
  VectorXd u1 = one_hot_vector(K, 1);
  VectorXd v1(K);
  v1 << 1, 0, 0, 0;

  EXPECT_EQ(K, u1.size());
  for (int i = 0; i < K; i++) {
    EXPECT_EQ(u1[i], v1[i]);
  }

  int k = 4;
  VectorXd u2 = one_hot_vector(K, k);
  VectorXd v2(K);
  v2 << 0, 0, 0, 1;

  EXPECT_EQ(k, u2.size());
  for (int i = 0; i < k; i++) {
    EXPECT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, one_hot_vector_throw) {
  using stan::math::one_hot_vector;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  int K = 5;
  int k = 2;

  EXPECT_THROW(one_hot_vector(K, K + 1), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, 0), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, -1), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, inf), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, nan), std::domain_error);

  EXPECT_THROW(one_hot_vector(-1, k), std::domain_error);
  EXPECT_THROW(one_hot_vector(inf, k), std::domain_error);
  EXPECT_THROW(one_hot_vector(nan, k), std::domain_error);
}
