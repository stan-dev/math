#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <gtest/gtest.h>
#include <CL/cl.hpp>
#include <algorithm>
#include <vector>

TEST(MathMatrixOpenCL, matrix_cl_copy_cache) {
  auto m = stan::math::matrix_d::Random(100, 100).eval();
  auto m_back = stan::math::matrix_d(100, 100);
  stan::math::matrix_cl d11(100, 100);
  stan::math::copy(d11, m);
  stan::math::copy(m_back, d11);
}

TEST(MathMatrixOpenCL, matrix_cl_var_copy) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> m(5, 5);
  double pos_ = 1.1;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      m(i, j) = pos_++;

  stan::math::matrix_cl d1_gpu(5, 5);
  stan::math::copy(d1_gpu, m);
  stan::math::matrix_d d1_cpu_return(5, 5);
  stan::math::copy(d1_cpu_return, d1_gpu);
  EXPECT_EQ(1.1, d1_cpu_return(0, 0));
  EXPECT_EQ(6.1, d1_cpu_return(1, 0));
  EXPECT_EQ(11.1, d1_cpu_return(2, 0));
  EXPECT_EQ(16.1, d1_cpu_return(3, 0));
  EXPECT_EQ(21.1, d1_cpu_return(4, 0));
}

TEST(MathMatrixOpenCL, matrix_cl_copy) {
  stan::math::vector_d d1;
  stan::math::vector_d d1_a;
  stan::math::vector_d d1_b;
  stan::math::matrix_d d2;
  stan::math::matrix_d d2_a;
  stan::math::matrix_d d2_b;
  stan::math::matrix_d d0;
  d1.resize(3);
  d1_a.resize(3);
  d1_b.resize(3);
  d1 << 1, 2, 3;
  d2.resize(2, 3);
  d2_a.resize(2, 3);
  d2_b.resize(2, 3);
  d2 << 1, 2, 3, 4, 5, 6;
  // vector
  stan::math::matrix_cl d11_cl(3, 1);
  stan::math::matrix_cl d111_cl(3, 1);
  EXPECT_NO_THROW(d11_cl = stan::math::to_matrix_cl(d1_cpu));
  EXPECT_NO_THROW(d111_cl = stan::math::copy_cl(d11_cl));
  EXPECT_NO_THROW(d1_a_cpu = stan::math::from_matrix_cl(d11_cl));
  EXPECT_NO_THROW(d1_b_cpu = stan::math::from_matrix_cl(d111_cl));
  EXPECT_EQ(1, d1_a_cpu(0));
  EXPECT_EQ(2, d1_a_cpu(1));
  EXPECT_EQ(3, d1_a_cpu(2));
  EXPECT_EQ(1, d1_b_cpu(0));
  EXPECT_EQ(2, d1_b_cpu(1));
  EXPECT_EQ(3, d1_b_cpu(2));
}

TEST(MathMatrixCL, matrix_cl_matrix_copy) {
  stan::math::matrix_d d2_cpu;
  stan::math::matrix_d d2_a_cpu;
  stan::math::matrix_d d2_b_cpu;
  stan::math::matrix_d d0_cpu;
  d2_cpu.resize(2, 3);
  d2_a_cpu.resize(2, 3);
  d2_b_cpu.resize(2, 3);
  d2_cpu << 1, 2, 3, 4, 5, 6;
  // matrix
  stan::math::matrix_cl d00_cl;
  stan::math::matrix_cl d000_cl;
  stan::math::matrix_cl d22_cl(2, 3);
  stan::math::matrix_cl d222_cl(2, 3);
  EXPECT_NO_THROW(d22_cl = stan::math::to_matrix_cl(d2_cpu));
  EXPECT_NO_THROW(d222_cl = stan::math::copy_cl(d22_cl));
  EXPECT_NO_THROW(d2_a_cpu = stan::math::from_matrix_cl(d22_cl));
  EXPECT_NO_THROW(d2_b_cpu = stan::math::from_matrix_cl(d222_cl));
  EXPECT_EQ(1, d2_a_cpu(0, 0));
  EXPECT_EQ(2, d2_a_cpu(0, 1));
  EXPECT_EQ(3, d2_a_cpu(0, 2));
  EXPECT_EQ(4, d2_a_cpu(1, 0));
  EXPECT_EQ(5, d2_a_cpu(1, 1));
  EXPECT_EQ(6, d2_a_cpu(1, 2));
  EXPECT_EQ(1, d2_b_cpu(0, 0));
  EXPECT_EQ(2, d2_b_cpu(0, 1));
  EXPECT_EQ(3, d2_b_cpu(0, 2));
  EXPECT_EQ(4, d2_b_cpu(1, 0));
  EXPECT_EQ(5, d2_b_cpu(1, 1));
  EXPECT_EQ(6, d2_b_cpu(1, 2));
  // zero sized copy
  EXPECT_NO_THROW(d00_cl = stan::math::to_matrix_cl(d0_cpu));
  EXPECT_NO_THROW(d0_cpu = stan::math::from_matrix_cl(d00_cl));
  EXPECT_NO_THROW(d000_cl = stan::math::copy_cl(d00_cl));
}

TEST(MathMatrixCL, matrix_cl_matrix_copy_arithmetic) {
  int test_val = 5;
  // Use this for successful copy
  stan::math::matrix_cl d22_cl(1, 1);
  // Use this for failure copy
  stan::math::matrix_cl d222_cl(2, 3);
  EXPECT_NO_THROW(d22_cl = stan::math::to_matrix_cl(test_val));
  EXPECT_NO_THROW(test_val = stan::math::from_matrix_cl<int>(d22_cl));
  EXPECT_THROW(d222_cl = stan::math::to_matrix_cl(test_val),
               std::invalid_argument);
  EXPECT_THROW(test_val = stan::math::from_matrix_cl<int>(d222_cl),
               std::invalid_argument);
}

TEST(MathMatrixCL, matrix_cl_pack_unpack_copy_lower) {
  int size = 42;
  int packed_size = size * (size + 1) / 2;
  std::vector<double> packed_mat(packed_size);
  std::vector<double> packed_mat_dst(packed_size);
  for (size_t i = 0; i < packed_mat.size(); i++) {
    packed_mat[i] = i;
  }
  stan::math::matrix_d m_flat_cpu(size, size);
  auto m_cl = stan::math::packed_copy<stan::math::TriangularViewCL::Lower>(
      packed_mat, size);
  m_flat_cpu = stan::math::from_matrix_cl(m_cl);
  size_t pos = 0;
  for (size_t j = 0; j < size; ++j) {
    for (size_t i = 0; i < j; i++) {
      EXPECT_EQ(m_flat_cpu(i, j), 0.0);
    }
    for (size_t i = j; i < size; ++i) {
      EXPECT_EQ(m_flat_cpu(i, j), packed_mat[pos]);
      pos++;
    }
  }
  packed_mat_dst
      = stan::math::packed_copy<stan::math::TriangularViewCL::Lower>(m_cl);
  for (size_t i = 0; i < packed_mat.size(); i++) {
    EXPECT_EQ(packed_mat[i], packed_mat_dst[i]);
  }
}

TEST(MathMatrixCL, matrix_cl_pack_unpack_copy_upper) {
  int size = 51;
  int packed_size = size * (size + 1) / 2;
  std::vector<double> packed_mat(packed_size);
  std::vector<double> packed_mat_dst(packed_size);
  for (size_t i = 0; i < packed_mat.size(); i++) {
    packed_mat[i] = i;
  }
  stan::math::matrix_d m_flat_cpu(size, size);
  auto m_cl = stan::math::packed_copy<stan::math::TriangularViewCL::Upper>(
      packed_mat, size);
  m_flat_cpu = stan::math::from_matrix_cl(m_cl);
  size_t pos = 0;
  for (size_t j = 0; j < size; ++j) {
    for (size_t i = 0; i <= j; i++) {
      EXPECT_EQ(m_flat_cpu(i, j), packed_mat[pos]);
      pos++;
    }
    for (size_t i = j + 1; i < size; ++i) {
      EXPECT_EQ(m_flat_cpu(i, j), 0.0);
    }
  }
  packed_mat_dst
      = stan::math::packed_copy<stan::math::TriangularViewCL::Upper>(m_cl);
  for (size_t i = 0; i < packed_mat.size(); i++) {
    EXPECT_EQ(packed_mat[i], packed_mat_dst[i]);
  }
}

TEST(MathMatrixCL, matrix_cl_pack_unpack_copy_exception) {
  std::vector<double> packed_mat;
  stan::math::matrix_cl m_cl_zero;
  EXPECT_NO_THROW(stan::math::packed_copy<stan::math::TriangularViewCL::Upper>(
      packed_mat, 0));
  EXPECT_NO_THROW(
      stan::math::packed_copy<stan::math::TriangularViewCL::Upper>(m_cl_zero));
  EXPECT_THROW(stan::math::packed_copy<stan::math::TriangularViewCL::Upper>(
                   packed_mat, 1),
               std::invalid_argument);
}
#endif
