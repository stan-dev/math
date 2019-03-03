#ifdef STAN_OPENCL

#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/diagonal_multiply.hpp>
#include <stan/math/opencl/copy.hpp>
#include <gtest/gtest.h>
#include <algorithm>

using stan::math::matrix_cl;
using stan::math::matrix_d;

TEST(MathMatrix, multiply_m_exception_pass_diagonal_mul) {
  auto m0 = matrix_d::Ones(3, 2).eval();
  auto m1 = matrix_d::Ones(1, 3).eval();
  matrix_cl mm0(m0);
  matrix_cl mm1(m1);
  stan::math::diagonal_multiply(mm0, 1.0);
  EXPECT_NO_THROW(stan::math::diagonal_multiply(mm0, 1.0));
  EXPECT_NO_THROW(stan::math::diagonal_multiply(mm1, 1.0));
}

TEST(MathMatrix, diagonal_multiply_value_check) {
  matrix_d m1(1, 1);
  matrix_d m1_dst(1, 1);
  m1 << 2;
  matrix_cl mm1(m1);
  matrix_cl mm1_dst(mm1.rows(), mm1.cols());
  EXPECT_NO_THROW(mm1_dst = stan::math::diagonal_multiply(mm1, 3.0));
  stan::math::copy(m1_dst, mm1_dst);
  EXPECT_EQ(6.0, m1_dst(0, 0));

  matrix_d m2_dst(1, 3);
  matrix_d m2 = matrix_d::Ones(1, 3) * 2;
  matrix_cl mm2(m2);
  matrix_cl mm2_dst(mm2.rows(), mm2.cols());
  EXPECT_NO_THROW(mm2_dst = stan::math::diagonal_multiply(mm2, 3.0));
  stan::math::copy(m2_dst, mm2_dst);
  EXPECT_EQ(6.0, m2_dst(0, 0));
  EXPECT_EQ(2.0, m2_dst(0, 1));
  EXPECT_EQ(2.0, m2_dst(0, 2));

  matrix_d m3_dst(3, 3);
  matrix_d m3 = matrix_d::Ones(3, 3) * 2;
  matrix_cl mm3(m3);
  matrix_cl mm3_dst(mm3.rows(), mm3.cols());
  EXPECT_NO_THROW(mm3_dst = stan::math::diagonal_multiply(mm3, 3.0));
  stan::math::copy(m3_dst, mm3_dst);
  EXPECT_EQ(6.0, m3_dst(0, 0));
  EXPECT_EQ(2.0, m3_dst(0, 1));
  EXPECT_EQ(2.0, m3_dst(0, 2));
  EXPECT_EQ(2.0, m3_dst(1, 0));
  EXPECT_EQ(6.0, m3_dst(1, 1));
  EXPECT_EQ(2.0, m3_dst(1, 2));
  EXPECT_EQ(2.0, m3_dst(2, 0));
  EXPECT_EQ(2.0, m3_dst(2, 1));
  EXPECT_EQ(6.0, m3_dst(2, 2));
}

#endif
