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
  stan::math::matrix_cl v_cl(v);
  stan::math::matrix_cl rv_cl(rv);
  stan::math::matrix_cl m_cl(1, 1);
  EXPECT_NO_THROW(m_cl = rv_cl * v_cl);
}

TEST(MathMatrix, one_dim_zero_matrix) {
  stan::math::matrix_d m0(5, 0);
  stan::math::matrix_d m1(0, 3);

  stan::math::matrix_cl m0_cl(m0);
  stan::math::matrix_cl m1_cl(m1);
  EXPECT_NO_THROW(m0_cl * m1_cl);

  EXPECT_NO_THROW(m0_cl * 2.0);
  EXPECT_NO_THROW(2.0 * m0_cl);

  EXPECT_NO_THROW(m1_cl * 2.0);
  EXPECT_NO_THROW(2.0 * m1_cl);
}

TEST(MathMatrix, zero_result_matrix) {
  stan::math::matrix_d m0(0, 5);
  stan::math::matrix_d m1(5, 0);

  stan::math::matrix_cl m0_cl(m0);
  stan::math::matrix_cl m1_cl(m1);
  EXPECT_NO_THROW(m0_cl * m1_cl);
}

TEST(MathMatrix, zero_size_input_matrix) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(0, 0);

  stan::math::matrix_cl m0_cl(m0);
  stan::math::matrix_cl m1_cl(m1);
  EXPECT_NO_THROW(m0_cl * m1_cl);

  EXPECT_NO_THROW(m0_cl * 2.0);
  EXPECT_NO_THROW(2.0 * m0_cl);
}

TEST(MathMatrix, non_matching_dim_excpetion) {
  stan::math::matrix_d m0(5, 3);
  stan::math::matrix_d m1(2, 6);

  stan::math::matrix_cl m0_cl(m0);
  stan::math::matrix_cl m1_cl(m1);
  EXPECT_THROW(m0_cl * m1_cl, std::invalid_argument);
}

TEST(MathMatrix, multiply_scalar) {
  auto v = stan::math::vector_d::Random(25).eval();
  stan::math::vector_d v_cl_res(25);
  auto rv = stan::math::row_vector_d::Random(25).eval();
  stan::math::row_vector_d rv_cl_res(25);
  auto m = stan::math::matrix_d::Random(5, 5).eval();
  stan::math::matrix_d m_cl_res(5, 5);

  stan::math::matrix_cl v_cl(v);
  v_cl = v_cl * 2.0;
  stan::math::copy(v_cl_res, v_cl);

  stan::math::matrix_cl rv_cl(rv);
  rv_cl = rv_cl * 2.0;
  stan::math::copy(rv_cl_res, rv_cl);

  stan::math::matrix_cl m_cl(m);
  m_cl = m_cl * 2.0;
  stan::math::copy(m_cl_res, m_cl);

  v = v * 2.0;
  rv = rv * 2.0;
  m = m * 2.0;

  EXPECT_MATRIX_NEAR(v, v_cl_res, 1e-10);
  EXPECT_MATRIX_NEAR(rv, rv_cl_res, 1e-10);
  EXPECT_MATRIX_NEAR(m, m_cl_res, 1e-10);
}

TEST(MathMatrix, row_vector_vector) {
  auto v = stan::math::vector_d::Random(5).eval();
  auto rv = stan::math::row_vector_d::Random(5).eval();
  stan::math::matrix_d m0(1, 1);
  stan::math::matrix_d m0_cl_res(1, 1);
  stan::math::matrix_d m1(5, 5);
  stan::math::matrix_d m1_cl_res(5, 5);

  m0 = rv * v;
  m1 = v * rv;

  stan::math::matrix_cl v_cl(v);
  stan::math::matrix_cl rv_cl(rv);
  stan::math::matrix_cl m0_cl(1, 1);
  stan::math::matrix_cl m1_cl(5, 5);

  m0_cl = rv_cl * v_cl;
  m1_cl = v_cl * rv_cl;
  stan::math::copy(m0_cl_res, m0_cl);
  stan::math::copy(m1_cl_res, m1_cl);

  EXPECT_MATRIX_NEAR(m0, m0_cl_res, 1e-10);
  EXPECT_MATRIX_NEAR(m1, m1_cl_res, 1e-10);
}

TEST(AgradRevMatrix, multiply_small) {
  using stan::math::multiply;
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  auto m2 = stan::math::matrix_d::Random(3, 3).eval();
  stan::math::matrix_d m3_cl_res(3, 3);

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = m11 * m22;

  stan::math::copy(m3_cl_res, m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}

TEST(AgradRevMatrix, multiply_big) {
  using stan::math::multiply;
  int size = 512;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  auto m2 = stan::math::matrix_d::Random(size, size).eval();
  stan::math::matrix_d m3_cl_res(size, size);

  stan::math::matrix_cl m11(m1);
  stan::math::matrix_cl m22(m2);

  auto m3 = (m1 * m2).eval();

  auto m33 = m11 * m22;

  stan::math::copy(m3_cl_res, m33);

  EXPECT_MATRIX_NEAR(m3, m3_cl_res, 1e-10);
}
#endif
#endif
