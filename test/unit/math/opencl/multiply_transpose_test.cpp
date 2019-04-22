#ifdef STAN_OPENCL

#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/multiply_transpose.hpp>
#include <stan/math/opencl/copy.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>
boost::random::mt19937 rng;

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrix, multiply_transpose_exception_fail_zero) {
  stan::math::row_vector_d rv(0);
  stan::math::matrix_d m(0, 3);
  stan::math::matrix_d m1(2, 0);
  stan::math::matrix_cl mm(m);
  stan::math::matrix_cl mm2(m1);
  stan::math::matrix_cl rvv(rv);
  stan::math::matrix_cl ans_mm1(1, 1);
  stan::math::matrix_cl ans_mm2(0, 0);
  stan::math::matrix_cl ans_mm3(2, 2);
  EXPECT_NO_THROW(ans_mm1 = stan::math::multiply_transpose(rvv));
  EXPECT_NO_THROW(ans_mm2 = stan::math::multiply_transpose(mm));
  EXPECT_NO_THROW(ans_mm3 = stan::math::multiply_transpose(mm2));
}

TEST(MathMatrix, multiply_m_m_exception_pass_dim) {
  stan::math::matrix_d m1(1, 3);
  stan::math::matrix_d m2(3, 5);
  stan::math::matrix_cl mm1(m1);
  stan::math::matrix_cl mm2(m2);
  stan::math::matrix_cl mm3a(1, 1);
  stan::math::matrix_cl mm3b(3, 3);
  EXPECT_NO_THROW(mm3a = stan::math::multiply_transpose(mm1));
  EXPECT_NO_THROW(mm3b = stan::math::multiply_transpose(mm2));
}

TEST(MathMatrix, multiply_zero_size) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;
  stan::math::matrix_cl v00(v0);
  stan::math::matrix_cl rv00(rv0);
  stan::math::matrix_cl m00(m0);
  EXPECT_NO_THROW(stan::math::multiply_transpose(v00));
  EXPECT_NO_THROW(stan::math::multiply_transpose(rv00));
  EXPECT_NO_THROW(stan::math::multiply_transpose(m00));
}

TEST(AgradRevMatrix, multiply_transposed_small) {
  using stan::math::multiply;
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m2(3, 3);
  stan::math::matrix_d m2_gpu(3, 3);

  m2 = m1 * m1.transpose();

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(3, 3);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);

  EXPECT_MATRIX_NEAR(m2, m2_gpu, 1e-10);
}

TEST(AgradRevMatrix, multiply_transposed_big) {
  using stan::math::multiply;
  stan::math::matrix_d m2, m2_gpu;

  int size = 500;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  m2.resize(size, size);
  m2_gpu.resize(size, size);
  m2 = m1 * m1.transpose();

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(size, size);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);
  EXPECT_MATRIX_NEAR(m2, m2_gpu, 1e-10);
}

TEST(AgradRevMatrix, multiply_transposed_big_non_square) {
  using stan::math::multiply;
  stan::math::matrix_d m2, m2_gpu;

  int size_x = 200;
  int size_y = 500;
  auto m1 = stan::math::matrix_d::Random(size_x, size_y).eval();
  m2.resize(size_x, size_x);
  m2_gpu.resize(size_x, size_x);

  m2 = m1 * m1.transpose();

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(size_x, size_x);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);

  EXPECT_MATRIX_NEAR(m2, m2_gpu, 1e-10);
}

TEST(AgradRevMatrix, multiply_transposed_big_x) {
  using stan::math::multiply;
  stan::math::matrix_d m2, m2_gpu;

  int size_x = 200;
  int size_y = 50;
  auto m1 = stan::math::matrix_d::Random(size_x, size_y).eval();
  m2.resize(size_x, size_x);
  m2_gpu.resize(size_x, size_x);

  m2 = m1 * m1.transpose();

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(size_x, size_x);
  m22 = stan::math::multiply_transpose(m11);
  stan::math::copy(m2_gpu, m22);

  EXPECT_MATRIX_NEAR(m2, m2_gpu, 1e-10);
}

#endif
