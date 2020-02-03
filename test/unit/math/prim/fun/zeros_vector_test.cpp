#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, zeros_vector) {
  for (int K = 0; K < 5; K++) {
    Eigen::VectorXd v = Eigen::VectorXd::Zero(K);
    expect_matrix_eq(v, stan::math::zeros_vector(K));
  }
}

TEST(MathFunctions, zeros_vector_throw) {
  EXPECT_THROW(stan::math::zeros_vector(-1), std::domain_error);
}
