#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, ones_row_vector) {
  for (int K = 0; K < 5; K++) {
    Eigen::RowVectorXd y = Eigen::RowVectorXd::Constant(K, 1);
    expect_matrix_eq(y, stan::math::ones_row_vector(K));
  }
}

TEST(MathFunctions, ones_row_vector_throw) {
  EXPECT_THROW(stan::math::ones_row_vector(-1), std::domain_error);
}
