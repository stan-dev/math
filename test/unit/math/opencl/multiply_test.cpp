#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#ifdef STAN_OPENCL
boost::random::mt19937 rng;

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrix, vector_row_vector) {
  stan::math::vector_d v(3);
  stan::math::row_vector_d rv(3);
  stan::math::matrix_cl v_gpu(v);
  stan::math::matrix_cl rv_gpu(rv);
  stan::math::matrix_cl m_gpu(1, 1);
  EXPECT_NO_THROW(m_gpu = rv_gpu * v_gpu);
}

TEST(MathMatrix, one_dim_zero_matrix) {
  stan::math::matrix_d m0(5, 0);
  stan::math::matrix_d m1(0, 3);

  stan::math::matrix_cl m0_gpu(m0);
  stan::math::matrix_cl m1_gpu(m1);
  EXPECT_NO_THROW(m0_gpu * m1_gpu);

  EXPECT_NO_THROW(m0_gpu * 2.0);
  EXPECT_NO_THROW(2.0 * m0_gpu);

  EXPECT_NO_THROW(m1_gpu * 2.0);
  EXPECT_NO_THROW(2.0 * m1_gpu);
}

TEST(MathMatrix, zero_result_matrix) {
  stan::math::matrix_d m0(0, 5);
  stan::math::matrix_d m1(5, 0);

  stan::math::matrix_cl m0_gpu(m0);
  stan::math::matrix_cl m1_gpu(m1);
  EXPECT_NO_THROW(m0_gpu * m1_gpu);
}

TEST(MathMatrix, zero_size_input_matrix) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(0, 0);

  stan::math::matrix_cl m0_gpu(m0);
  stan::math::matrix_cl m1_gpu(m1);
  EXPECT_NO_THROW(m0_gpu * m1_gpu);

  EXPECT_NO_THROW(m0_gpu * 2.0);
  EXPECT_NO_THROW(2.0 * m0_gpu);
}

TEST(MathMatrix, non_matching_dim_excpetion) {
  stan::math::matrix_d m0(5, 3);
  stan::math::matrix_d m1(2, 6);

  stan::math::matrix_cl m0_gpu(m0);
  stan::math::matrix_cl m1_gpu(m1);
  EXPECT_THROW(m0_gpu * m1_gpu, std::invalid_argument);
}

TEST(MathMatrix, multiply_scalar) {
  auto v = stan::math::vector_d::Random(25).eval();
  stan::math::vector_d v_gpu_res(25);
  auto rv = stan::math::row_vector_d::Random(25).eval();
  stan::math::row_vector_d rv_gpu_res(25);
  auto m = stan::math::matrix_d::Random(5, 5).eval();
  stan::math::matrix_d m_gpu_res(5, 5);

  stan::math::matrix_cl v_gpu(v);
  v_gpu = v_gpu * 2.0;
  stan::math::copy(v_gpu_res, v_gpu);

  stan::math::matrix_cl rv_gpu(rv);
  rv_gpu = rv_gpu * 2.0;
  stan::math::copy(rv_gpu_res, rv_gpu);

  stan::math::matrix_cl m_gpu(m);
  m_gpu = m_gpu * 2.0;
  stan::math::copy(m_gpu_res, m_gpu);

  v = v * 2.0;
  rv = rv * 2.0;
  m = m * 2.0;

  EXPECT_MATRIX_NEAR(v, v_gpu_res, 1e-10);
  EXPECT_MATRIX_NEAR(rv, rv_gpu_res, 1e-10);
  EXPECT_MATRIX_NEAR(m, m_gpu_res, 1e-10);
}

TEST(MathMatrix, row_vector_vector) {
  auto v = stan::math::vector_d::Random(5).eval();
  auto rv = stan::math::row_vector_d::Random(5).eval();
  stan::math::matrix_d m0(1, 1);
  stan::math::matrix_d m0_gpu_res(1, 1);
  stan::math::matrix_d m1(5, 5);
  stan::math::matrix_d m1_gpu_res(5, 5);

  m0 = rv * v;
  m1 = v * rv;

  stan::math::matrix_cl v_gpu(v);
  stan::math::matrix_cl rv_gpu(rv);
  stan::math::matrix_cl m0_gpu(1, 1);
  stan::math::matrix_cl m1_gpu(5, 5);

  m0_gpu = rv_gpu * v_gpu;
  m1_gpu = v_gpu * rv_gpu;
  stan::math::copy(m0_gpu_res, m0_gpu);
  stan::math::copy(m1_gpu_res, m1_gpu);

  EXPECT_MATRIX_NEAR(m0, m0_gpu_res, 1e-10);
  EXPECT_MATRIX_NEAR(m1, m1_gpu_res, 1e-10);
}

TEST(AgradRevMatrix, multiply_small) {
  using stan::math::multiply;
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_gpu_res(3, 3);

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = m11 * m22;

  stan::math::copy(m3_gpu_res, m33);

  EXPECT_MATRIX_NEAR(m3, m3_gpu_res, 1e-10);
}

TEST(AgradRevMatrix, multiply_big) {
  using stan::math::multiply;
  int size = 512;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_gpu_res(size, size);

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = m11 * m22;

  stan::math::copy(m3_gpu_res, m33);

  EXPECT_MATRIX_NEAR(m3, m3_gpu_res, 1e-10);
}
#endif
#endif
