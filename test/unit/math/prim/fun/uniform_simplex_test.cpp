#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, uniform_simplex) {
  using Eigen::VectorXd;
  using stan::math::uniform_simplex;

  VectorXd u0 = uniform_simplex(0);
  EXPECT_EQ(0, u0.size());

  VectorXd u1 = uniform_simplex(1);
  VectorXd v1(1);
  v1 << 1;
  EXPECT_EQ(1, u1.size());
  EXPECT_EQ(u1[0], v1[0]);

  int K = 4;
  VectorXd u2 = uniform_simplex(K);
  VectorXd v2(K);
  v2 << 0.25, 0.25, 0.25, 0.25;
  EXPECT_EQ(K, u2.size());
  for (int i = 0; i < K; i++) {
    EXPECT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, uniform_simplex_throw) {
  using stan::math::uniform_simplex;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(uniform_simplex(-1), std::domain_error);
  EXPECT_THROW(uniform_simplex(inf), std::domain_error);
  EXPECT_THROW(uniform_simplex(nan), std::domain_error);
}
