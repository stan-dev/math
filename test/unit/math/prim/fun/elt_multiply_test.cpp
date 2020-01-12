#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, eltMultiplyVec) {
  stan::math::vector_d v1(2);
  stan::math::vector_d v2(2);
  v1 << 1, 2;
  v2 << 10, 100;
  stan::math::vector_d v = stan::math::elt_multiply(v1, v2);
  EXPECT_FLOAT_EQ(10.0, v(0));
  EXPECT_FLOAT_EQ(200.0, v(1));
}
TEST(MathMatrixPrimMat, eltMultiplyVecException) {
  stan::math::vector_d v1(2);
  stan::math::vector_d v2(3);
  v1 << 1, 2;
  v2 << 10, 100, 1000;
  EXPECT_THROW(stan::math::elt_multiply(v1, v2), std::invalid_argument);
}
TEST(MathMatrixPrimMat, eltMultiplyRowVec) {
  stan::math::row_vector_d v1(2);
  stan::math::row_vector_d v2(2);
  v1 << 1, 2;
  v2 << 10, 100;
  stan::math::row_vector_d v = stan::math::elt_multiply(v1, v2);
  EXPECT_FLOAT_EQ(10.0, v(0));
  EXPECT_FLOAT_EQ(200.0, v(1));
}
TEST(MathMatrixPrimMat, eltMultiplyRowVecException) {
  stan::math::row_vector_d v1(2);
  stan::math::row_vector_d v2(3);
  v1 << 1, 2;
  v2 << 10, 100, 1000;
  EXPECT_THROW(stan::math::elt_multiply(v1, v2), std::invalid_argument);
}
TEST(MathMatrixPrimMat, eltMultiplyMatrix) {
  stan::math::matrix_d m1(2, 3);
  stan::math::matrix_d m2(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  m2 << 10, 100, 1000, 10000, 100000, 1000000;
  stan::math::matrix_d m = stan::math::elt_multiply(m1, m2);

  EXPECT_EQ(2, m.rows());
  EXPECT_EQ(3, m.cols());
  EXPECT_FLOAT_EQ(10.0, m(0, 0));
  EXPECT_FLOAT_EQ(200.0, m(0, 1));
  EXPECT_FLOAT_EQ(3000.0, m(0, 2));
  EXPECT_FLOAT_EQ(40000.0, m(1, 0));
  EXPECT_FLOAT_EQ(500000.0, m(1, 1));
  EXPECT_FLOAT_EQ(6000000.0, m(1, 2));
}
TEST(MathMatrixPrimMat, eltMultiplyMatrixException) {
  stan::math::matrix_d m1(2, 3);
  stan::math::matrix_d m2(2, 4);
  stan::math::matrix_d m3(4, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  m2 << -1, -2, -3, -4, -5, -6, -7, -8;
  m3 << 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24;
  EXPECT_THROW(stan::math::elt_multiply(m1, m2), std::invalid_argument);
  EXPECT_THROW(stan::math::elt_multiply(m1, m3), std::invalid_argument);
}
