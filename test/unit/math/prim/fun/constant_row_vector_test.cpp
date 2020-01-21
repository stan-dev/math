#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, constant_row_vector) {
  using Eigen::RowVectorXd;
  using stan::math::constant_row_vector;

  RowVectorXd u0 = constant_row_vector(0, 1);
  EXPECT_EQ(0, u0.size());

  double c = 3.14;
  for (int K = 1; K < 5; K++) {
    RowVectorXd v = RowVectorXd::Constant(K, c);
    expect_matrix_eq(v, constant_row_vector(K, c));
  }
}

TEST(MathFunctions, constant_row_vector_throw) {
  EXPECT_THROW(stan::math::constant_row_vector(-1, 3.14), std::domain_error);
}
