#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <CL/cl2.hpp>
#include <algorithm>
#include <vector>

TEST(MathMatrixGPU, copy_rvalue_constructor) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  int N = 100;
  MatrixXd a = MatrixXd::Random(N, N);
  MatrixXd b = MatrixXd::Random(N, N);
  matrix_cl<double> c_cl(a + b);  // the problem
  // attempt to scramble the memory that was used by the temporary
  MatrixXd w = a - b;

  // make compiler think a and b were changed, so second sum can not be
  // optimized away
  volatile int i = 19;
  volatile double no_opt = a(i);
  a(i) = no_opt;
  no_opt = b(i);
  b(i) = no_opt;

  MatrixXd correct = a + b;
  MatrixXd result = from_matrix_cl(c_cl);
  for (int i = 0; i < N * N; i++) {
    EXPECT_EQ(result(i), correct(i));
  }
  no_opt = w.array().sum();
}

TEST(MathMatrixGPU, copy_rvalue_constructor_eval) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  int N = 100;
  MatrixXd a = MatrixXd::Random(N, N);
  MatrixXd b = MatrixXd::Random(N, N);
  matrix_cl<double> c_cl(
      (a + b).eval());  // the problem, this case with .eval()
  // attempt to scramble the memory that was used by the temporary
  MatrixXd w = a - b;

  // make compiler think a and b were changed, so second sum can not be
  // optimized away
  volatile int i = 19;
  volatile double no_opt = a(i);
  a(i) = no_opt;
  no_opt = b(i);
  b(i) = no_opt;

  MatrixXd correct = a + b;
  MatrixXd result = from_matrix_cl(c_cl);
  for (int i = 0; i < N * N; i++) {
    EXPECT_EQ(result(i), correct(i));
  }
  no_opt = w.array().sum();
}

TEST(MathMatrixGPU, copy_rvalue_to_matrix_cl) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  int N = 100;
  MatrixXd a = MatrixXd::Random(N, N);
  MatrixXd b = MatrixXd::Random(N, N);
  matrix_cl<double> c_cl = stan::math::to_matrix_cl(a + b);  // the problem
  // attempt to scramble the memory that was used by the temporary
  MatrixXd w = a - b;

  // make compiler think a and b were changed, so second sum can not be
  // optimized away
  volatile int i = 19;
  volatile double no_opt = a(i);
  a(i) = no_opt;
  no_opt = b(i);
  b(i) = no_opt;

  MatrixXd correct = a + b;
  MatrixXd result = from_matrix_cl(c_cl);
  for (int i = 0; i < N * N; i++) {
    EXPECT_EQ(result(i), correct(i));
  }
  no_opt = w.array().sum();
}

TEST(MathMatrixGPU, copy_rvalue_to_matrix_cl_eval) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  int N = 100;
  MatrixXd a = MatrixXd::Random(N, N);
  MatrixXd b = MatrixXd::Random(N, N);
  matrix_cl<double> c_cl
      = stan::math::to_matrix_cl((a + b).eval());  // the problem
  // attempt to scramble the memory that was used by the temporary
  MatrixXd w = a - b;

  // make compiler think a and b were changed, so second sum can not be
  // optimized away
  volatile int i = 19;
  volatile double no_opt = a(i);
  a(i) = no_opt;
  no_opt = b(i);
  b(i) = no_opt;

  MatrixXd correct = a + b;
  MatrixXd result = from_matrix_cl(c_cl);
  for (int i = 0; i < N * N; i++) {
    EXPECT_EQ(result(i), correct(i));
  }
  no_opt = w.array().sum();
}

TEST(MathMatrixGPU, copy_rvalue_to_matrix_cl_scalar) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  double a = 1.0;
  double b = 2.0;
  double c = 3.0;
  double d = 4.0;
  matrix_cl<double> c_cl
      = stan::math::to_matrix_cl(a + b + c + d);  // the problem
  // attempt to scramble the memory that was used by the temporary
  double f = a - b - c - d;

  // make compiler think a and b were changed, so second sum can not be
  // optimized away
  volatile double no_opt = a;
  a = no_opt;
  no_opt = b;
  b = no_opt;

  double correct = a + b + c + d;
  MatrixXd result = from_matrix_cl(c_cl);
  EXPECT_EQ(result(0), correct);
  no_opt = f + b + c + d;
}

TEST(MathMatrixGPU, copy_rvalue_constructor_scalar) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  double a = 1.0;
  double b = 2.0;
  double c = 3.0;
  double d = 4.0;
  matrix_cl<double> c_cl(a + b + c + d);  // the problem
  // attempt to scramble the memory that was used by the temporary
  double f = a - b - c - d;

  // make compiler think a and b were changed, so second sum can not be
  // optimized away
  volatile double no_opt = a;
  a = no_opt;
  no_opt = b;
  b = no_opt;

  double correct = a + b + c + d;
  MatrixXd result = from_matrix_cl(c_cl);
  EXPECT_EQ(result(0), correct);
  no_opt = f + b + c + d;
}

TEST(MathMatrixGPU, std_vector_copy_rvalue_constructor) {
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;

  matrix_cl<double> c_cl(
      std::vector<double>({1, 2, 3, 4, 5, 6}));  // the problem

  std::vector<int> correct{1, 2, 3, 4, 5, 6};
  MatrixXd result = from_matrix_cl(c_cl);
  for (int i = 0; i < correct.size(); i++) {
    EXPECT_EQ(result(i), correct[i]);
  }
}

TEST(MathMatrixGPU, std_vector_copy_rvalue_to_matrix_cl) {
  using Eigen::MatrixXd;
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;

  matrix_cl<double> c_cl = stan::math::to_matrix_cl(
      std::vector<double>({1, 2, 3, 4, 5, 6}));  // the problem

  std::vector<int> correct{1, 2, 3, 4, 5, 6};
  MatrixXd result = from_matrix_cl(c_cl);
  for (int i = 0; i < correct.size(); i++) {
    EXPECT_EQ(result(i), correct[i]);
  }
}

TEST(MathMatrixGPU, matrix_cl_vector_copy) {
  stan::math::vector_d d1_cpu;
  stan::math::vector_d d1_a_cpu;
  stan::math::vector_d d1_b_cpu;
  d1_cpu.resize(3);
  d1_a_cpu.resize(3);
  d1_b_cpu.resize(3);
  d1_cpu << 1, 2, 3;
  // vector
  stan::math::matrix_cl<double> d11_cl(3, 1);
  stan::math::matrix_cl<double> d111_cl(3, 1);
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

TEST(MathMatrixGPU, matrix_cl_std_vector_copy) {
  std::vector<double> d1_cpu{10.0, 20.0, 30.0};
  std::vector<double> d_empty;
  std::vector<double> d1_cpu_ret;

  stan::math::matrix_cl<double> d11_cl = stan::math::to_matrix_cl(d1_cpu);
  stan::math::matrix_cl<double> d_empty_cl;
  EXPECT_NO_THROW(d_empty_cl = stan::math::to_matrix_cl(d_empty));
  EXPECT_NO_THROW(d1_cpu_ret
                  = stan::math::from_matrix_cl<std::vector<double>>(d11_cl));
  EXPECT_STD_VECTOR_EQ(d1_cpu_ret, d1_cpu);
}

TEST(MathMatrixGPU, matrix_cl_std_vector_eigen_copy) {
  Eigen::VectorXd v1(3);
  v1 << 1, 2, 3;
  Eigen::VectorXd v2(3);
  v2 << 4, 5, 6;
  std::vector<Eigen::VectorXd> v{v1, v2};
  stan::math::matrix_cl<double> v_cl;
  EXPECT_NO_THROW(v_cl = stan::math::to_matrix_cl(v));
  std::vector<Eigen::VectorXd> v_ret;
  EXPECT_NO_THROW(
      v_ret = stan::math::from_matrix_cl<std::vector<Eigen::VectorXd>>(v_cl));
  EXPECT_MATRIX_EQ(v_ret[0], v1);
  EXPECT_MATRIX_EQ(v_ret[1], v2);
}

TEST(MathMatrixGPU, matrix_cl_packed_std_vector_copy_rvalue) {
  stan::math::matrix_d d1_cpu_ret;

  stan::math::matrix_cl<double> d11_cl
      = stan::math::packed_copy<stan::math::matrix_cl_view::Lower>(
          std::vector<double>({10.0, 20.0, 30.0}), 2);
  EXPECT_NO_THROW(d1_cpu_ret = stan::math::from_matrix_cl(d11_cl));
  EXPECT_EQ(10.0, d1_cpu_ret(0, 0));
  EXPECT_EQ(20.0, d1_cpu_ret(1, 0));
  EXPECT_EQ(30.0, d1_cpu_ret(1, 1));
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
  stan::math::matrix_cl<double> d00_cl;
  stan::math::matrix_cl<double> d000_cl;
  stan::math::matrix_cl<double> d22_cl(2, 3);
  stan::math::matrix_cl<double> d222_cl(2, 3);
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
  double test_val = 5;
  // Use this for successful copy
  stan::math::matrix_cl<double> d22_cl(1, 1);
  EXPECT_NO_THROW(d22_cl = stan::math::to_matrix_cl(test_val));
  EXPECT_NO_THROW(test_val = stan::math::from_matrix_cl<double>(d22_cl));
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
  stan::math::matrix_cl<double> m_cl
      = stan::math::packed_copy<stan::math::matrix_cl_view::Lower>(packed_mat,
                                                                   size);
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
  packed_mat_dst = stan::math::packed_copy(m_cl);
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
  auto m_cl = stan::math::packed_copy<stan::math::matrix_cl_view::Upper>(
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
  packed_mat_dst = stan::math::packed_copy(m_cl);
  for (size_t i = 0; i < packed_mat.size(); i++) {
    EXPECT_EQ(packed_mat[i], packed_mat_dst[i]);
  }
}

TEST(MathMatrixCL, matrix_cl_pack_unpack_copy_exception) {
  std::vector<double> packed_mat;
  stan::math::matrix_cl<double> m_cl_zero;
  EXPECT_NO_THROW(stan::math::packed_copy<stan::math::matrix_cl_view::Upper>(
      packed_mat, 0));
  m_cl_zero.view(stan::math::matrix_cl_view::Upper);
  EXPECT_NO_THROW(stan::math::packed_copy(m_cl_zero));
  m_cl_zero.view(stan::math::matrix_cl_view::Entire);
  EXPECT_THROW(stan::math::packed_copy(m_cl_zero), std::invalid_argument);
  EXPECT_THROW(
      stan::math::packed_copy<stan::math::matrix_cl_view::Upper>(packed_mat, 1),
      std::invalid_argument);
}
#endif
