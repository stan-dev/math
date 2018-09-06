#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/multiply.hpp>
#include <stan/math/gpu/copy.hpp>
#include <stan/math/gpu/lower_tri_inverse.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixGPU, inverse_gpu_exception) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_gpu m2(m1);
  stan::math::matrix_gpu m3(2, 3);
  using stan::math::lower_triangular_inverse;
  EXPECT_THROW(m3 = lower_triangular_inverse(m2), std::invalid_argument);
}

TEST(MathMatrixGPU, inverse_gpu_small) {
  auto m1 = stan::math::matrix_d::Random(3, 3).eval();
  m1.triangularView<Eigen::StrictlyUpper>() = stan::math::matrix_d::Zero(3, 3).eval();

  stan::math::matrix_d m1_cpu(3, 3);
  stan::math::matrix_d m1_cl(3, 3);

  m1_cpu = stan::math::mdivide_left_tri<Eigen::Lower>(m1);

  stan::math::matrix_gpu m2(m1);
  stan::math::matrix_gpu m3(3, 3);

  EXPECT_NO_THROW(m3 = stan::math::lower_triangular_inverse(m2));
  stan::math::copy(m1_cl, m3);

  EXPECT_MATRIX_NEAR(m1_cl, m1_cpu, 1e-8);
}

TEST(MathMatrixGPU, inverse_gpu_big) {
  int size = 128;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  m1.triangularView<Eigen::StrictlyUpper>() = stan::math::matrix_d::Zero(size, size).eval();
  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);

  m1_cpu = stan::math::mdivide_left_tri<Eigen::Lower>(m1);

  stan::math::matrix_gpu m2(m1);
  stan::math::matrix_gpu m3(size, size);

  EXPECT_NO_THROW(m3 = stan::math::lower_triangular_inverse(m2));
  EXPECT_NO_THROW(stan::math::copy(m1_cl, m3));

  EXPECT_MATRIX_NEAR(m1_cl, m1_cpu, 1e-8);

  /*
  for (int i = 0; i < m1.rows(); i++) {
    for (int j = 0; j < m1.cols(); j++) {
      std::cout << m1(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;std::cout << std::endl;
  for (int i = 0; i < m1_cpu.rows(); i++) {
    for (int j = 0; j < m1_cpu.cols(); j++) {
      std::cout << m1_cpu(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;std::cout << std::endl;
  for (int i = 0; i < m2_cpu.rows(); i++) {
    for (int j = 0; j < m2_cpu.cols(); j++) {
      std::cout << m2_cpu(i, j) << " ";
      EXPECT_NEAR(m1_cpu(i, j), m2_cpu(i, j), 1e-8);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;std::cout << std::endl;
  */
  //EXPECT_MATRIX_NEAR(m2_cpu, m1_cpu, 1e-8);
}
#endif
