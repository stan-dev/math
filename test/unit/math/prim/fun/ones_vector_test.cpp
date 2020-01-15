#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, ones_vector) {
  using Eigen::VectorXd;
  using stan::math::ones_vector;

  VectorXd u0 = ones_vector(0);
  EXPECT_EQ(0, u0.size());

  VectorXd u1 = ones_vector(1);
  VectorXd v1(1);
  v1 << 1;
  EXPECT_EQ(1, u1.size());
  EXPECT_EQ(u1[0], v1[0]);

  int K = 4;
  VectorXd u2 = ones_vector(K);
  VectorXd v2(K);
  v2 << 1, 1, 1, 1;
  EXPECT_EQ(K, u2.size());
  for (int i = 0; i < K; i++) {
    EXPECT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, ones_vector_throw) {
  using stan::math::ones_vector;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(ones_vector(-1), std::domain_error);
  EXPECT_THROW(ones_vector(inf), std::domain_error);
  EXPECT_THROW(ones_vector(nan), std::domain_error);
}
