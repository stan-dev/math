#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/basic_matrix_gpu.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, copy_triangular_m_exception_pass) {
  stan::math::matrix_gpu m0;

  EXPECT_NO_THROW(stan::math::copy_triangular(m0, stan::math::UPPER));
  EXPECT_NO_THROW(stan::math::copy_triangular(m0, stan::math::LOWER));

  stan::math::matrix_gpu m1(1, 1);

  EXPECT_NO_THROW(stan::math::copy_triangular(m1, stan::math::UPPER));
  EXPECT_NO_THROW(stan::math::copy_triangular(m1, stan::math::LOWER));
}

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

TEST(MathMatrixGPU, copy_triangular_m_pass) {
  stan::math::matrix_d m0(2, 2);
  stan::math::matrix_d m0_dst(2, 2);
  m0 << 1, 2,
        3, 4;
  m0_dst << 0, 0,
            0, 0;

  stan::math::matrix_gpu m00(m0);
  stan::math::matrix_gpu m00_dst(m0_dst);

  EXPECT_NO_THROW(m00_dst = stan::math::copy_triangular(m00,
    stan::math::UPPER));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m00_dst));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(2, m0_dst(0, 1));
  EXPECT_EQ(0, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));

  EXPECT_NO_THROW(m00_dst = stan::math::copy_triangular(m00,
    stan::math::LOWER));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m00_dst));
  EXPECT_EQ(1, m0_dst(0, 0));
  EXPECT_EQ(0, m0_dst(0, 1));
  EXPECT_EQ(3, m0_dst(1, 0));
  EXPECT_EQ(4, m0_dst(1, 1));
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

TEST(MathMatrixGPU, copy_submatrix_exception) {
  stan::math::matrix_d d1;

  d1.resize(3, 3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(2, 2);
  EXPECT_NO_THROW(stan::math::copy_submatrix(d11,
    d22, 1, 1, 0, 0, 2, 2));


  EXPECT_THROW(stan::math::copy_submatrix(d11, d22,
    1, 1, 0, 0, 4, 4), std::domain_error);
  EXPECT_THROW(stan::math::copy_submatrix(d11, d22,
    4, 4, 0, 0, 2, 2), std::domain_error);
  EXPECT_THROW(stan::math::copy_submatrix(d11, d22,
    1, 1, 3, 3, 4, 4), std::domain_error);
  EXPECT_THROW(stan::math::copy_submatrix(d11, d22,
    1, 1, 3, 3, 3, 3), std::domain_error);
}
#endif
