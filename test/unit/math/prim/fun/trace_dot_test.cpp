#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, trace_dot) {
  using stan::math::matrix_d;
  using stan::math::trace_dot;

  matrix_d a(2, 3);
  matrix_d b(3, 2);
  a << 1, 2, 3, 4, 5, 6;
  b << 7, 8, 9, 10, 11, 12;

  // trace(A * B) = trace([[58,64],[139,154]]) = 58 + 154 = 212
  EXPECT_FLOAT_EQ(212, trace_dot(a, b));
}

TEST(MathMatrixPrim, trace_dot_square) {
  using stan::math::matrix_d;
  using stan::math::trace_dot;

  matrix_d a(2, 2);
  matrix_d b(2, 2);
  a << 1, 2, 3, 4;
  b << 5, 6, 7, 8;

  // trace(A * B) = trace([[19,22],[43,50]]) = 19 + 50 = 69
  EXPECT_FLOAT_EQ(69, trace_dot(a, b));
}

TEST(MathMatrixPrim, trace_dot_1x1) {
  using stan::math::matrix_d;
  using stan::math::trace_dot;

  matrix_d a(1, 1);
  matrix_d b(1, 1);
  a << 3;
  b << 7;

  EXPECT_FLOAT_EQ(21, trace_dot(a, b));
}

TEST(MathMatrixPrim, trace_dot_size_zero) {
  using stan::math::matrix_d;
  using stan::math::trace_dot;

  matrix_d a00, b00;
  EXPECT_FLOAT_EQ(0, trace_dot(a00, b00));
}

TEST(MathMatrixPrim, trace_dot_dimension_mismatch) {
  using stan::math::matrix_d;
  using stan::math::trace_dot;

  matrix_d a(2, 3);
  matrix_d b(2, 3);
  a << 1, 2, 3, 4, 5, 6;
  b << 1, 2, 3, 4, 5, 6;

  EXPECT_THROW(trace_dot(a, b), std::invalid_argument);
}
