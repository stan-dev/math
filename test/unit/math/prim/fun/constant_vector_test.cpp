#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, constant_vector) {
  using Eigen::VectorXd;
  using stan::math::constant_vector;

  VectorXd u0 = constant_vector(0, 1);
  EXPECT_EQ(0, u0.size());

  double c = 3.14;
  VectorXd u1 = constant_vector(1, c);
  VectorXd v1(1);
  v1 << c;
  EXPECT_EQ(1, u1.size());
  EXPECT_EQ(u1[0], v1[0]);

  int k = 4;
  VectorXd u2 = constant_vector(k, c);
  VectorXd v2(k);
  v2 << c, c, c, c;
  EXPECT_EQ(k, u2.size());
  for (int i = 0; i < k; i++) {
    EXPECT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, constant_vector_throw) {
  using stan::math::constant_vector;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double c = 3.14;

  EXPECT_THROW(constant_vector(-1, c), std::domain_error);
  EXPECT_THROW(constant_vector(inf, c), std::domain_error);
  EXPECT_THROW(constant_vector(nan, c), std::domain_error);
}
