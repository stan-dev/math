#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/copy.hpp>
#include <stan/math/gpu/cholesky_decompose.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrix, cholesky_decompose) {
  stan::math::matrix_d m0(3, 3);
  m0 <<  25, 15, -5, 15, 18,  0, -5,  0, 11;

  stan::math::matrix_d m1(4, 4);
  m1 <<  18, 22,  54,  42, 22, 70,  86,  62, 54, 86, 174, 134, 42, 62, 134, 106;
  
  stan::math::matrix_gpu m0_gpu(m0);
  stan::math::matrix_gpu m1_gpu(m1);
  
  stan::math::matrix_d m0_res = stan::math::cholesky_decompose(m0);
  stan::math::matrix_d m1_res = stan::math::cholesky_decompose(m1);

  stan::math::matrix_gpu m0_chol_gpu = stan::math::cholesky_decompose(m0_gpu);
  stan::math::matrix_gpu m1_chol_gpu = stan::math::cholesky_decompose(m1_gpu);
  
  stan::math::copy(m0, m0_chol_gpu);
  stan::math::copy(m1, m1_chol_gpu);

  EXPECT_MATRIX_NEAR(m0, m0_res, 1e-8);
  EXPECT_MATRIX_NEAR(m1, m1_res, 1e-8);  
}

void cholesky_decompose_test(int size) {
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_pos_def = m1 * m1.transpose() + size*Eigen::MatrixXd::Identity(size, size);

  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);

  m1_cpu = stan::math::cholesky_decompose(m1_pos_def);

  stan::math::matrix_gpu m2(m1_pos_def);
  auto m3 = stan::math::cholesky_decompose(m2);
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

TEST(MathMatrix, cholesky_decompose_small) {
  cholesky_decompose_test(10);
  cholesky_decompose_test(50);
  cholesky_decompose_test(100);
}

TEST(MathMatrix, cholesky_decompose_big) {
  cholesky_decompose_test(500);
  cholesky_decompose_test(1000);
  cholesky_decompose_test(2000);
}


#endif
