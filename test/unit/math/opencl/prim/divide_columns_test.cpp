#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/prim/divide_columns.hpp>
#include <gtest/gtest.h>
#include <algorithm>

using stan::math::divide_columns;
using stan::math::from_matrix_cl;
using stan::math::matrix_cl;
using stan::math::matrix_d;
using stan::math::row_vector_d;
using stan::math::vector_d;

TEST(MathMatrixCL, divide_columns_v_exception_pass) {
  vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_NO_THROW(divide_columns(d11, d22));
}

TEST(MathMatrixCL, divide_columns_v_exception_pass_zero) {
  vector_d d1, d2;
  d1.resize(0);
  d2.resize(0);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_NO_THROW(divide_columns(d11, d22));
}

TEST(MathMatrixCL, divide_columns_v_exception_pass_invalid_arg) {
  row_vector_d d1, d2;

  d1.resize(2);
  d2.resize(3);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_THROW(divide_columns(d11, d22), std::invalid_argument);
}

TEST(MathMatrixCL, divide_columns_rv_exception_pass) {
  row_vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_NO_THROW(divide_columns(d11, d22));
}

TEST(MathMatrixCL, divide_columns_rv_exception_pass_zero) {
  row_vector_d d1, d2;

  d1.resize(0);
  d2.resize(0);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_NO_THROW(divide_columns(d11, d22));
}

TEST(MathMatrixCL, divide_columns_rv_exception_fail_invalid_arg) {
  row_vector_d d1, d2;

  d1.resize(2);
  d2.resize(3);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_THROW(divide_columns(d11, d22), std::invalid_argument);
}

TEST(MathMatrixCL, divide_columns_m_exception_pass_simple) {
  matrix_d d1, d2;

  d1.resize(2, 3);
  d2.resize(2, 3);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_THROW(divide_columns(d11, d22), std::invalid_argument);
}

TEST(MathMatrixCL, divide_columns_m_exception_pass_zero) {
  matrix_d d1, d2;
  d1.resize(0, 0);
  d2.resize(0, 0);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_NO_THROW(divide_columns(d11, d22));
}

TEST(MathMatrixCL, divide_columns_m_exception_fail_invalid_arg) {
  matrix_d d1, d2;
  d1.resize(2, 3);
  d2.resize(3, 3);
  matrix_cl<double> d11(d1);
  matrix_cl<double> d22(d2);
  EXPECT_THROW(divide_columns(d11, d22), std::invalid_argument);
}

TEST(MathMatrixCL, divide_columns_non_matching_sizes_exception) {
  vector_d v1(2);
  v1 << 1, 2;
  vector_d v2(3);
  v2 << 10, 100, 1000;

  row_vector_d rv1(2);
  rv1 << 1, 2;
  row_vector_d rv2(3);
  rv2 << 10, 100, 1000;

  matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  matrix_d m2(3, 2);
  m2 << 10, 100, 1000, 0, -10, -12;

  matrix_cl<double> v11(v1);
  matrix_cl<double> v22(v2);
  matrix_cl<double> rv11(rv1);
  matrix_cl<double> rv22(rv2);
  matrix_cl<double> m11(m1);
  matrix_cl<double> m22(m2);

  EXPECT_THROW(divide_columns(v11, v22), std::invalid_argument);
  EXPECT_THROW(divide_columns(rv11, rv22), std::invalid_argument);
  EXPECT_THROW(divide_columns(m11, m22), std::invalid_argument);
}

TEST(MathMatrixCL, divide_columns_value_vector_check) {
  vector_d v1(3);
  vector_d v2(3);
  v1 << 10, 20, 30;
  v2 << 10, 20, 30;
  vector_d v3(3);

  matrix_cl<double> v11(v1);
  matrix_cl<double> v22(v2);
  matrix_cl<double> v33(3, 1);

  EXPECT_NO_THROW(divide_columns(v11, v22));
  v3 = from_matrix_cl(v11);
  EXPECT_EQ(1, v3(0));
  EXPECT_EQ(1, v3(1));
  EXPECT_EQ(1, v3(2));
}

TEST(MathMatrixCL, divide_columns_value_vector_d_check) {
  row_vector_d rv1(3);
  row_vector_d rv2(3);
  rv1 << 10, 20, 30;
  rv2 << 10, 20, 30;
  row_vector_d rv3(3);

  matrix_cl<double> rv11(rv1);
  matrix_cl<double> rv22(rv2);
  matrix_cl<double> rv33(1, 3);

  EXPECT_NO_THROW(divide_columns(rv11, rv22));
  rv3 = from_matrix_cl(rv11);
  EXPECT_EQ(1, rv3(0));
  EXPECT_EQ(1, rv3(1));
  EXPECT_EQ(1, rv3(2));
}

TEST(MathMatrixCL, divide_columns_value_matrix_check) {
  matrix_d m1(3, 3);
  matrix_d m2(3, 1);
  m1 << 10, 20, 30, 20, 40, 60, 30, 60, 90;
  m2 << 10, 10, 10;
  matrix_d m3(3, 3);

  matrix_cl<double> m11(m1);
  matrix_cl<double> m22(m2);
  matrix_cl<double> m33(3, 3);

  EXPECT_NO_THROW(divide_columns(m11, m22));
  m3 = from_matrix_cl(m11);
  EXPECT_EQ(1, m3(0, 0));
  EXPECT_EQ(2, m3(0, 1));
  EXPECT_EQ(3, m3(0, 2));
  EXPECT_EQ(2, m3(1, 0));
  EXPECT_EQ(4, m3(1, 1));
  EXPECT_EQ(6, m3(1, 2));
  EXPECT_EQ(3, m3(2, 0));
  EXPECT_EQ(6, m3(2, 1));
  EXPECT_EQ(9, m3(2, 2));
}

#endif
