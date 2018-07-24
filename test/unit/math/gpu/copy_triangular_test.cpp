#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/copy_triangular.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, copy_triangular_m_exception_pass) {
  stan::math::matrix_gpu m0;

  EXPECT_NO_THROW(stan::math::copy_triangular<stan::math::gpu::Upper>(m0));
  EXPECT_NO_THROW(stan::math::copy_triangular<stan::math::gpu::Lower>(m0));

  stan::math::matrix_gpu m1(1, 1);

  EXPECT_NO_THROW(stan::math::copy_triangular<stan::math::gpu::Upper>(m1));
  EXPECT_NO_THROW(stan::math::copy_triangular<stan::math::gpu::Lower>(m1));
}

TEST(MathMatrixGPU, copy_triangular_m_pass) {
  stan::math::matrix_d m0(2, 2);
  stan::math::matrix_d m0_dst(2, 2);
  m0 << 1, 2, 3, 4;
  m0_dst << 0, 0, 0, 0;

  stan::math::matrix_gpu m00(m0);
  stan::math::matrix_gpu m00_dst(m0_dst);

  EXPECT_NO_THROW(m00_dst
                  = stan::math::copy_triangular<stan::math::gpu::Upper>(m00));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m00_dst));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(0, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));

  EXPECT_NO_THROW(m00_dst
                  = stan::math::copy_triangular<stan::math::gpu::Lower>(m00));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m00_dst));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(0, m0_dst(0, 1));
  EXPECT_EQ(3, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));
}
#endif
