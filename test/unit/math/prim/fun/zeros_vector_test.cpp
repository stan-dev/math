#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, zeros_vector) {
  using Eigen::VectorXd;
  using stan::math::zeros_vector;

  VectorXd u0 = zeros_vector(0);
  EXPECT_EQ(0, u0.size());

  VectorXd u1 = zeros_vector(1);
  VectorXd v1(1);
  v1 << 0;
  EXPECT_EQ(1, u1.size());
  EXPECT_EQ(u1[0], v1[0]);

  int K = 4;
  VectorXd u2 = zeros_vector(K);
  VectorXd v2(K);
  v2 << 0, 0, 0, 0;
  EXPECT_EQ(K, u2.size());
  for (int i = 0; i < K; i++) {
    EXPECT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, zeros_vector_throw) {
  using stan::math::zeros_vector;

  EXPECT_THROW(zeros_vector(-1), std::domain_error);
}
