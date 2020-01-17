#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, zeros_vector) {
  using Eigen::VectorXd;
  using stan::math::zeros_vector;

  VectorXd u0 = zeros_vector(0);
  EXPECT_EQ(0, u0.size());

  for (int K = 1; K < 5; K++) {
    Eigen::VectorXd v = Eigen::VectorXd::Zero(K);
    expect_matrix_eq(v, zeros_vector(K));
  }
}

TEST(MathFunctions, zeros_vector_throw) {
  EXPECT_THROW(stan::math::zeros_vector(-1), std::domain_error);
}
