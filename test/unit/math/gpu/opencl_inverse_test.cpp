#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/basic_matrix_gpu.hpp>
#include <stan/math/gpu/multiply_matrix_gpu.hpp>
#include <stan/math/gpu/inverse_gpu.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, inverse_gpu_exception) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_gpu m2(m1);
  stan::math::matrix_gpu m3(2, 3);
  using stan::math::lower_triangular_inverse;
  EXPECT_THROW(m3 = lower_triangular_inverse(m2), std::invalid_argument);
}

TEST(MathMatrixGPU, inverse_gpu_small) {
  stan::math::matrix_d m1(3, 3);
  stan::math::matrix_d m1a(3, 3);
  stan::math::matrix_d m1_dst(3, 3);
  m1 << 1, 0, 0, 2, 3, 0, 4, 5, 6;

  m1a = stan::math::inverse(m1);

  stan::math::matrix_gpu m2(m1);
  stan::math::matrix_gpu m3(3, 3);

  EXPECT_NO_THROW(m3 = stan::math::lower_triangular_inverse(m2));
  EXPECT_NO_THROW(stan::math::copy(m1_dst, m3));
  for (int i = 0; i < 3; i++) {
  for (int j = 0; j < 3; j++)
    if (i == j)
      EXPECT_NEAR(m1_dst(i, j), m1a(i, j), 1e-10);
    else
      EXPECT_NEAR(m1_dst(i, j), m1a(i, j), 1e-10);
  }
}

TEST(MathMatrixGPU, inverse_gpu_big) {
  int size = 512;
  stan::math::matrix_d m1(size, size);
  stan::math::matrix_d m1a(size, size);
  stan::math::matrix_d m1_dst(size, size);

  for (int i = 0; i < size; i++)
  for (int j = 0; j < size; j++) {
    if ( j <= i )
      m1(i, j) = (j%10)+1;
    else
      m1(i, j) = 0.0;
  }
  m1a = stan::math::inverse(m1);

  stan::math::matrix_gpu m2(m1);
  stan::math::matrix_gpu m3(size, size);

  EXPECT_NO_THROW(m3 = stan::math::lower_triangular_inverse(m2));
  EXPECT_NO_THROW(stan::math::copy(m1_dst, m3));

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      EXPECT_NEAR(m1_dst(i, j), m1a(i, j), 1e-10);
  }
}
#endif
