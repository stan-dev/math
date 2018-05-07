#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, matrix_gpu_creation) {
  stan::math::vector_d d1;
  stan::math::matrix_d d2;

  d1.resize(3);
  d2.resize(2, 3);
  EXPECT_NO_THROW(stan::math::matrix_gpu A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_gpu d11(d1));
  EXPECT_NO_THROW(stan::math::matrix_gpu d22(d2));
}

TEST(MathMatrixGPU, matrix_gpu_copy) {
  stan::math::vector_d d1;
  stan::math::vector_d d1_a;
  stan::math::vector_d d1_b;
  stan::math::matrix_d d2;
  stan::math::matrix_d d2_a;
  stan::math::matrix_d d2_b;

  d1.resize(3);
  d1_a.resize(3);
  d1_b.resize(3);
  d1 << 1, 2, 3;
  d2.resize(2, 3);
  d2_a.resize(2, 3);
  d2_b.resize(2, 3);

  d2 << 1, 2, 3,
        4, 5, 6;
  // vector
  stan::math::matrix_gpu d11(3, 1);
  stan::math::matrix_gpu d111(3, 1);
  EXPECT_NO_THROW(stan::math::copy(d1, d11));
  EXPECT_NO_THROW(stan::math::copy(d11, d111));
  EXPECT_NO_THROW(stan::math::copy(d11, d1_a));
  EXPECT_NO_THROW(stan::math::copy(d111, d1_b));
  EXPECT_EQ(1, d1_a(0));
  EXPECT_EQ(2, d1_a(1));
  EXPECT_EQ(3, d1_a(2));
  EXPECT_EQ(1, d1_b(0));
  EXPECT_EQ(2, d1_b(1));
  EXPECT_EQ(3, d1_b(2));
  // matrix
  stan::math::matrix_gpu d22(2, 3);
  stan::math::matrix_gpu d222(2, 3);
  EXPECT_NO_THROW(stan::math::copy(d2, d22));
  EXPECT_NO_THROW(stan::math::copy(d22, d222));
  EXPECT_NO_THROW(stan::math::copy(d22, d2_a));
  EXPECT_NO_THROW(stan::math::copy(d222, d2_b));
  EXPECT_EQ(1, d2_a(0, 0));
  EXPECT_EQ(2, d2_a(0, 1));
  EXPECT_EQ(3, d2_a(0, 2));
  EXPECT_EQ(4, d2_a(1, 0));
  EXPECT_EQ(5, d2_a(1, 1));
  EXPECT_EQ(6, d2_a(1, 2));
  EXPECT_EQ(1, d2_b(0, 0));
  EXPECT_EQ(2, d2_b(0, 1));
  EXPECT_EQ(3, d2_b(0, 2));
  EXPECT_EQ(4, d2_b(1, 0));
  EXPECT_EQ(5, d2_b(1, 1));
  EXPECT_EQ(6, d2_b(1, 2));
}
