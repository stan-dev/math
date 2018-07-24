#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/diagonal_multiply.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#ifdef STAN_OPENCL
TEST(MathMatrix, multiply_m_exception_pass_diagonal_mul) {
  stan::math::matrix_d m0;
  m0.resize(3, 2);
  m0 << 1, 1, 1,
  1, 1, 1;
  stan::math::matrix_d m1;
  m1.resize(1, 3);
  m1 << 1, 1, 1;
  stan::math::matrix_gpu mm0(m0);
  stan::math::matrix_gpu mm1(m1);
  stan::math::diagonal_multiply(mm0, 1.0);
  EXPECT_NO_THROW(stan::math::diagonal_multiply(mm0, 1.0));
  EXPECT_NO_THROW(stan::math::diagonal_multiply(mm1, 1.0));
}

TEST(MathMatrix, diagonal_multiply_value_check) {
  stan::math::matrix_d m1;
  stan::math::matrix_d m1_dst;
  m1.resize(1, 1);
  m1_dst.resize(1, 1);
  m1 << 2;
  stan::math::matrix_gpu mm1(m1);
  stan::math::matrix_gpu mm1_dst(mm1.rows(),mm1.cols());
  EXPECT_NO_THROW(mm1_dst = stan::math::diagonal_multiply(mm1, 3.0));
  stan::math::copy(m1_dst, mm1_dst);
  EXPECT_EQ(6.0, m1_dst(0,0));
  
  stan::math::matrix_d m2;
  stan::math::matrix_d m2_dst;
  m2.resize(1, 3);
  m2_dst.resize(1, 3);
  m2 << 2, 2, 2;
  stan::math::matrix_gpu mm2(m2);
  stan::math::matrix_gpu mm2_dst(mm2.rows(),mm2.cols());
  EXPECT_NO_THROW(mm2_dst = stan::math::diagonal_multiply(mm2, 3.0));
  stan::math::copy(m2_dst, mm2_dst);
  EXPECT_EQ(6.0, m2_dst(0,0));
  EXPECT_EQ(2.0, m2_dst(0,1));
  EXPECT_EQ(2.0, m2_dst(0,2));
  
  stan::math::matrix_d m3;
  stan::math::matrix_d m3_dst;
  m3.resize(3, 3);
  m3_dst.resize(3, 3);
  m3 << 2, 2, 2, 2, 2, 2, 2, 2, 2;
  stan::math::matrix_gpu mm3(m3);
  stan::math::matrix_gpu mm3_dst(mm3.rows(),mm3.cols());
  EXPECT_NO_THROW(mm3_dst = stan::math::diagonal_multiply(mm3, 3.0));
  stan::math::copy(m3_dst, mm3_dst);
  EXPECT_EQ(6.0, m3_dst(0,0));
  EXPECT_EQ(2.0, m3_dst(0,1));
  EXPECT_EQ(2.0, m3_dst(0,2));
  EXPECT_EQ(2.0, m3_dst(1,0));
  EXPECT_EQ(6.0, m3_dst(1,1));
  EXPECT_EQ(2.0, m3_dst(1,2));
  EXPECT_EQ(2.0, m3_dst(2,0));
  EXPECT_EQ(2.0, m3_dst(2,1));
  EXPECT_EQ(6.0, m3_dst(2,2));
}

#endif
