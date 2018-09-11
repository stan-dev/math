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
/*
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
}*/

TEST(MathMatrixGPU, inverse_gpu_big) {
  int size = 64*4*4*8;
  boost::random::mt19937 rng;
  std::cout << "Size: " << size << std::endl;
  auto m1 = stan::math::matrix_d::Random(size, size).eval();
  m1.triangularView<Eigen::StrictlyUpper>() = stan::math::matrix_d::Zero(size, size).eval();
  for(int i = 0; i < m1.rows(); i++){
    for(int j = 0; j <= i;j++){
      m1(i, j) = stan::math::uniform_rng(-20000, 20000, rng);
    }
  }
  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);
  
  auto start_cpu = std::chrono::steady_clock::now();
  m1_cpu = stan::math::mdivide_left_tri<Eigen::Lower>(m1).eval();
  auto stop_cpu = std::chrono::steady_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop_cpu - start_cpu).count()
            << "ms" << std::endl;

  auto start_gpu = std::chrono::steady_clock::now();
  stan::math::matrix_gpu m3(size, size);
  stan::math::matrix_gpu m2(m1);  
  EXPECT_NO_THROW(m3 = stan::math::lower_triangular_inverse(m2));
  EXPECT_NO_THROW(stan::math::copy(m1_cl, m3));
  auto stop_gpu = std::chrono::steady_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop_gpu - start_gpu).count()
            << "ms" << std::endl;
  /*for(int i = 0; i < m1.rows();i++){
    for(int j = 0; j < m1.cols();j++){
      std::cout << m1_cpu(i,j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/
  /*
  for(int i = 0; i < m1.rows();i++){
    for(int j = 0; j < m1.cols();j++){
      std::cout << m1_cl(i,j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/
  double max_error=0;
  for(int i = 0; i < m1.rows(); i++){
    for(int j = 0; j < m1.cols();j++){
      double abs_err = abs(m1_cpu(i, j) - m1_cl(i, j));
      double a = std::max(abs_err/m1_cpu(i, j), abs_err/m1_cl(i, j));
      if( a > 1e-6){
        std::cout << i << ", " << j << a << std::endl;
      }
      max_error = std::max(max_error, a);
      
    }    
  }
  std::cout << "Max error: " << max_error << std::endl;
  //EXPECT_MATRIX_NEAR(m1_cpu, m1_cl, 1e-4); 
    
}
#endif
