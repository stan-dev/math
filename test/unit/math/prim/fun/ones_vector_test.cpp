#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, ones_vector) {
  using Eigen::VectorXd;
  using stan::math::ones_vector;

  VectorXd u0 = ones_vector(0);
  EXPECT_EQ(0, u0.size());

  for (int K = 1; K < 5; K++) {
    Eigen::VectorXd y = Eigen::VectorXd::Constant(K, 1);
    expect_matrix_eq(y, ones_vector(K));
  }
}

TEST(MathFunctions, ones_vector_throw) {
  EXPECT_THROW(stan::math::ones_vector(-1), std::domain_error);
}
