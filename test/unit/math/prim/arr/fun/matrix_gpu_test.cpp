#include <stan/math/prim/mat.hpp>
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
  stan::math::matrix_d d2;

  d1.resize(3);
  d2.resize(2, 3);
  stan::math::matrix_gpu d11(3, 1);
  stan::math::matrix_gpu d111(3, 1);
  stan::math::matrix_gpu d22(2, 3);
  EXPECT_NO_THROW(stan::math::copy(d1, d11));
  EXPECT_NO_THROW(stan::math::copy(d2, d22));
  EXPECT_NO_THROW(stan::math::copy(d11, d1));
  EXPECT_NO_THROW(stan::math::copy(d22, d2));
  EXPECT_NO_THROW(stan::math::copy(d11, d111));
}

