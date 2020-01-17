#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, constant_vector) {
  using Eigen::VectorXd;
  using stan::math::constant_vector;

  VectorXd u0 = constant_vector(0, 1);
  EXPECT_EQ(0, u0.size());

  double c = 3.14;
  for (int K = 1; K < 5; K++) {
    VectorXd v = VectorXd::Constant(K, c);
    expect_matrix_eq(v, constant_vector(K, c));
  }
}

TEST(MathFunctions, constant_vector_throw) {
  EXPECT_THROW(stan::math::constant_vector(-1, 3.14), std::domain_error);
}
