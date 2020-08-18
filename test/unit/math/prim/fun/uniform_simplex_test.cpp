#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, uniform_simplex) {
  using Eigen::VectorXd;
  using stan::math::uniform_simplex;

  for (int K = 1; K < 5; K++) {
    Eigen::VectorXd v = Eigen::VectorXd::Constant(K, 1.0 / K);
    EXPECT_MATRIX_FLOAT_EQ(v, uniform_simplex(K));
  }
}

TEST(MathFunctions, uniform_simplex_throw) {
  using stan::math::uniform_simplex;
  EXPECT_THROW(uniform_simplex(-1), std::domain_error);
  EXPECT_THROW(uniform_simplex(0), std::domain_error);
}
