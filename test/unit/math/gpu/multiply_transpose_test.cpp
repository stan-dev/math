#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/multiply_transpose.hpp>
#include <stan/math/gpu/copy.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#ifdef STAN_OPENCL
TEST(MathMatrix, multiply_transpose_exception_fail_zero) {
  stan::math::row_vector_d rv;
  stan::math::matrix_d m;
  stan::math::matrix_d m1;
  rv.resize(0);
  m.resize(0, 3);
  m1.resize(2, 0);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu mm2(m1);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu ans_mm1(1, 1);
  stan::math::matrix_gpu ans_mm2(0, 0);
  stan::math::matrix_gpu ans_mm3(2, 2);
  EXPECT_NO_THROW(ans_mm1 = stan::math::multiply_transpose(rvv));
  EXPECT_NO_THROW(ans_mm2 = stan::math::multiply_transpose(mm));
  EXPECT_NO_THROW(ans_mm3 = stan::math::multiply_transpose(mm2));
}

TEST(MathMatrix, multiply_m_m_exception_pass_dim) {
  stan::math::matrix_d m1, m2;

  m1.resize(1, 3);
  m2.resize(3, 5);
  stan::math::matrix_gpu mm1(m1);
  stan::math::matrix_gpu mm2(m2);
  stan::math::matrix_gpu mm3a(1, 1);
  stan::math::matrix_gpu mm3b(3, 3);
  EXPECT_NO_THROW(mm3a = stan::math::multiply_transpose(mm1));
  EXPECT_NO_THROW(mm3b = stan::math::multiply_transpose(mm2));
}

TEST(MathMatrix, multiply_zero_size) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  stan::math::matrix_gpu v00(v0);
  stan::math::matrix_gpu rv00(rv0);
  stan::math::matrix_gpu m00(m0);
  EXPECT_NO_THROW(stan::math::multiply_transpose(v00));
  EXPECT_NO_THROW(stan::math::multiply_transpose(rv00));
  EXPECT_NO_THROW(stan::math::multiply_transpose(m00));
}

TEST(AgradRevMatrix, multiply_transposed_small) {
  using stan::math::multiply;
  stan::math::matrix_d m1, m2, m2_gpu;

  m1.resize(3, 3);
  m2.resize(3, 3);
  m2_gpu.resize(3, 3);

  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  m2 = m1 * m1.transpose();

  stan::math::matrix_gpu m11(m1);
  stan::math::matrix_gpu m22(3, 3);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_NEAR(m2(i, j), m2_gpu(i, j), 1e-10);
}

TEST(AgradRevMatrix, multiply_transposed_big) {
  boost::random::mt19937 rng;
  using stan::math::multiply;
  stan::math::matrix_d m1, m2, m2_gpu;

  int size = 500;
  m1.resize(size, size);
  m2.resize(size, size);
  m2_gpu.resize(size, size);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) {
      m1(i, j) = stan::math::normal_rng(0.0, 1.0, rng);
    }
  m2 = m1 * m1.transpose();

  stan::math::matrix_gpu m11(m1);
  stan::math::matrix_gpu m22(size, size);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      EXPECT_NEAR(m2(i, j), m2_gpu(i, j), 1e-10);
}

TEST(AgradRevMatrix, multiply_transposed_big_non_square) {
  boost::random::mt19937 rng;
  using stan::math::multiply;
  stan::math::matrix_d m1, m2, m2_gpu;

  int size_x = 200;
  int size_y = 500;
  m1.resize(size_x, size_y);
  m2.resize(size_x, size_x);
  m2_gpu.resize(size_x, size_x);
  for (int i = 0; i < size_x; i++)
    for (int j = 0; j < size_y; j++) {
      m1(i, j) = stan::math::normal_rng(0.0, 1.0, rng);
    }

  m2 = m1 * m1.transpose();

  stan::math::matrix_gpu m11(m1);
  stan::math::matrix_gpu m22(size_x, size_x);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);

  for (int i = 0; i < size_x; i++)
    for (int j = 0; j < size_x; j++)
      EXPECT_NEAR(m2(i, j), m2_gpu(i, j), 1e-10);
}

#endif
