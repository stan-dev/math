#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, zeros_row_vector) {
  using Eigen::RowVectorXd;
  using stan::math::zeros_row_vector;

  RowVectorXd u0 = zeros_row_vector(0);
  EXPECT_EQ(0, u0.size());

  for (int K = 1; K < 5; K++) {
    Eigen::RowVectorXd v = Eigen::RowVectorXd::Zero(K);
    expect_matrix_eq(v, zeros_row_vector(K));
  }
}

TEST(MathFunctions, zeros_row_vector_throw) {
  EXPECT_THROW(stan::math::zeros_row_vector(-1), std::domain_error);
}
