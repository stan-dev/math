#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, minus) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  EXPECT_EQ(0, stan::math::minus(v0).size());
  EXPECT_EQ(0, stan::math::minus(rv0).size());
  EXPECT_EQ(0, stan::math::minus(m0).size());
}

TEST(MathFunctions, minus_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::minus(b));
}
