#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/add.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, add_v_exception_pass) {
  stan::math::vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(3, 1);
  EXPECT_NO_THROW(d33 = d11 + d22);
}

TEST(MathMatrixGPU, add_v_exception_pass_zero) {
  stan::math::vector_d d1, d2;
  d1.resize(0);
  d2.resize(0);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(0, 1);
  EXPECT_NO_THROW(d33 = d11 + d22);
}

TEST(MathMatrixGPU, add_v_exception_pass_invalid_arg) {
  stan::math::row_vector_d d1, d2;

  d1.resize(2);
  d2.resize(3);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(3, 0);
  EXPECT_THROW(d33 = d11 + d22, std::invalid_argument);
}

TEST(MathMatrixGPU, add_rv_exception_pass) {
  stan::math::row_vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(1, 3);
  EXPECT_NO_THROW(d33 = d11 + d22);
}

TEST(MathMatrixGPU, add_rv_exception_pass_zero) {
  stan::math::row_vector_d d1, d2;

  d1.resize(0);
  d2.resize(0);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(1, 0);
  EXPECT_NO_THROW(d33 = d11 + d22);
}

TEST(MathMatrixGPU, add_rv_exception_fail_invalid_arg) {
  stan::math::row_vector_d d1, d2;

  d1.resize(2);
  d2.resize(3);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(3, 1);
  EXPECT_THROW(d33 = d11 + d22, std::invalid_argument);
}

TEST(MathMatrixGPU, add_m_exception_pass_simple) {
  stan::math::matrix_d d1, d2;

  d1.resize(2, 3);
  d2.resize(2, 3);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(2, 3);
  EXPECT_NO_THROW(d33 = d11 + d22);
}

TEST(MathMatrixGPU, add_m_exception_pass_zero) {
  stan::math::matrix_d d1, d2;
  d1.resize(0, 0);
  d2.resize(0, 0);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(0, 0);
  EXPECT_NO_THROW(d33 = d11 + d22);
}

TEST(MathMatrixGPU, add_m_exception_fail_invalid_arg) {
  stan::math::matrix_d d1, d2;
  d1.resize(2, 3);
  d2.resize(3, 3);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  stan::math::matrix_cl d33(2, 3);
  EXPECT_THROW(d33 = d11 + d22, std::invalid_argument);
}

TEST(MathMatrixGPU, add_non_matching_sizes_exception) {
  stan::math::vector_d v1(2);
  v1 << 1, 2;
  stan::math::vector_d v2(3);
  v2 << 10, 100, 1000;

  stan::math::row_vector_d rv1(2);
  rv1 << 1, 2;
  stan::math::row_vector_d rv2(3);
  rv2 << 10, 100, 1000;

  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_d m2(3, 2);
  m2 << 10, 100, 1000, 0, -10, -12;

  using stan::math::add;
  using stan::math::matrix_cl;
  matrix_cl v11(v1);
  matrix_cl v22(v2);
  matrix_cl v33(v1);
  matrix_cl rv11(rv1);
  matrix_cl rv22(rv2);
  matrix_cl rv33(rv1);
  matrix_cl m11(m1);
  matrix_cl m22(m2);
  matrix_cl m33(m1);

  EXPECT_THROW(v33 = v11 + v22, std::invalid_argument);
  EXPECT_THROW(rv33 = rv11 + rv22, std::invalid_argument);
  EXPECT_THROW(m33 = m11 + m22, std::invalid_argument);
}

TEST(MathMatrixGPU, add_value_check) {
  stan::math::vector_d v1(3);
  v1 << 1, 2, 3;
  stan::math::vector_d v2(3);
  v2 << 10, 100, 1000;
  stan::math::vector_d v3(3);

  stan::math::row_vector_d rv1(3);
  rv1 << 1, 2, 3;
  stan::math::row_vector_d rv2(3);
  rv2 << 10, 100, 1000;
  stan::math::row_vector_d rv3(3);

  stan::math::matrix_d m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_d m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;
  stan::math::matrix_d m3(3, 3);

  using stan::math::add;
  using stan::math::matrix_cl;
  matrix_cl v11(v1);
  matrix_cl v22(v2);
  matrix_cl v33(3, 1);
  matrix_cl rv11(rv1);
  matrix_cl rv22(rv2);
  matrix_cl rv33(1, 3);
  matrix_cl m11(m1);
  matrix_cl m22(m2);
  matrix_cl m33(3, 3);

  EXPECT_NO_THROW(v33 = v11 + v22);
  EXPECT_NO_THROW(rv33 = rv11 + rv22);
  EXPECT_NO_THROW(m33 = m11 + m22);

  stan::math::copy(v3, v33);
  EXPECT_EQ(11, v3(0));
  EXPECT_EQ(102, v3(1));
  EXPECT_EQ(1003, v3(2));

  stan::math::copy(rv3, rv33);
  EXPECT_EQ(11, rv3(0));
  EXPECT_EQ(102, rv3(1));
  EXPECT_EQ(1003, rv3(2));

  stan::math::copy(m3, m33);
  EXPECT_EQ(11, m3(0, 0));
  EXPECT_EQ(102, m3(0, 1));
  EXPECT_EQ(1003, m3(0, 2));
  EXPECT_EQ(4, m3(1, 0));
  EXPECT_EQ(-5, m3(1, 1));
  EXPECT_EQ(-6, m3(1, 2));
  EXPECT_EQ(9, m3(2, 0));
  EXPECT_EQ(12, m3(2, 1));
  EXPECT_EQ(17, m3(2, 2));
}
#endif
