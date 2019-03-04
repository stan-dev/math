#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/lower_tri_inverse.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixGPU, inverse_gpu_exception) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl m2(m1);
  stan::math::matrix_cl m3(2, 3);
  using stan::math::lower_triangular_inverse;
  EXPECT_THROW(m3 = lower_triangular_inverse(m2), std::invalid_argument);
}

void inverse_test(int size) {
  boost::random::mt19937 rng;
  auto m1 = stan::math::matrix_d(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      m1(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    m1(i, i) = 10000.0;
    for (int j = i + 1; j < size; j++) {
      m1(i, j) = 0.0;
    }
  }
  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);

  m1_cpu = stan::math::mdivide_left_tri<Eigen::Lower>(m1);

  stan::math::matrix_cl m2(m1);
  auto m3 = stan::math::lower_triangular_inverse(m2);
  stan::math::copy(m1_cl, m3);
  double max_error = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j <= i; j++) {
      double abs_err = std::fabs(m1_cpu(i, j) - m1_cl(i, j));
      double a = std::max(abs_err / m1_cpu(i, j), abs_err / m1_cl(i, j));
      max_error = std::max(max_error, a);
    }
  }
  EXPECT_LT(max_error, 1e-8);
}
TEST(MathMatrixGPU, inverse_gpu_small) { inverse_test(3); }

TEST(MathMatrixGPU, inverse_gpu_1under_block_size) { inverse_test(31); }

TEST(MathMatrixGPU, inverse_gpu_1over_block_size) { inverse_test(33); }

TEST(MathMatrixGPU, inverse_gpu_big) { inverse_test(512); }

TEST(MathMatrixGPU, inverse_gpu_big_2) { inverse_test(700); }

TEST(MathMatrixGPU, inverse_gpu_big_3) { inverse_test(1500); }
#endif
