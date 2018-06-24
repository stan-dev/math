#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/copy_triangular_transposed_opencl.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, copy_triangular_transposed_m_exception_pass) {
  stan::math::matrix_gpu m1(1, 1);
  stan::math::matrix_gpu m0;

  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(m0,
    stan::math::LOWER_TO_UPPER_TRIANGULAR));
  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(m0,
    stan::math::LOWER_TO_UPPER_TRIANGULAR));
  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(m1,
    stan::math::LOWER_TO_UPPER_TRIANGULAR));
  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(m1,
    stan::math::UPPER_TO_LOWER_TRIANGULAR));
}

TEST(MathMatrixGPU, copy_triangular_transposed_m_pass) {
  stan::math::matrix_d m0(2, 2);
  stan::math::matrix_d m0_dst(2, 2);
  m0 << 1, 2,
        3, 4;

  stan::math::matrix_gpu m00(m0);
  stan::math::matrix_gpu m11(m0);

  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(m00,
    stan::math::LOWER_TO_UPPER_TRIANGULAR));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m00));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(3, m0_dst(0, 1));
  EXPECT_EQ(3, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));

  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(m11,
    stan::math::UPPER_TO_LOWER_TRIANGULAR));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m11));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(2, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));
}

#endif
