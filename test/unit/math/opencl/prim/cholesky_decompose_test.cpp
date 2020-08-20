#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/opencl.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixOpenCLPrim, cholesky_decompose_cl_exceptions) {
  using stan::math::matrix_cl;
  matrix_cl<double> m0;
  EXPECT_NO_THROW(stan::math::cholesky_decompose(m0));

  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  matrix_cl<double> m1_cl(m1);
  EXPECT_THROW(stan::math::cholesky_decompose(m1_cl), std::invalid_argument);

  // not pos-def
  stan::math::matrix_d m22(2, 2);
  m22 << 1.0, 2.0, 2.0, 3.0;
  matrix_cl<double> m22_cl(m22);
  EXPECT_THROW(stan::math::cholesky_decompose(m22_cl), std::domain_error);

  // not symmetric
  stan::math::matrix_d m_not_sym(2, 2);
  m_not_sym << 1.0, 2.0, 3.0, 4.0;
  matrix_cl<double> m_not_sym_cl(m_not_sym);
  EXPECT_THROW(stan::math::cholesky_decompose(m_not_sym_cl), std::domain_error);
}

TEST(MathMatrixOpenCLPrim, cholesky_decompose_non_inplace_cpu_vs_cl_small) {
  stan::math::matrix_d m0(3, 3);
  m0 << 25, 15, -5, 15, 18, 0, -5, 0, 11;

  stan::math::matrix_d m1(4, 4);
  m1 << 18, 22, 54, 42, 22, 70, 86, 62, 54, 86, 174, 134, 42, 62, 134, 106;

  stan::math::matrix_cl<double> m0_cl(m0);
  stan::math::matrix_cl<double> m1_cl(m1);

  stan::math::matrix_d m0_res = stan::math::cholesky_decompose(m0);
  stan::math::matrix_d m1_res = stan::math::cholesky_decompose(m1);

  stan::math::matrix_cl<double> m0_cl_res
      = stan::math::cholesky_decompose(m0_cl);
  stan::math::matrix_cl<double> m1_cl_res
      = stan::math::cholesky_decompose(m1_cl);

  m0 = stan::math::from_matrix_cl(m0_cl_res);
  m1 = stan::math::from_matrix_cl(m1_cl_res);

  EXPECT_MATRIX_NEAR(m0, m0_res, 1e-8);
  EXPECT_MATRIX_NEAR(m1, m1_res, 1e-8);
}

namespace {
void cholesky_decompose_test(int size) {
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_pos_def
      = m1 * m1.transpose() + size * Eigen::MatrixXd::Identity(size, size);

  stan::math::check_symmetric("cholesky_decompose", "m", m1_pos_def);
  Eigen::LLT<stan::math::matrix_d> llt(m1_pos_def.rows());
  llt.compute(m1_pos_def);
  stan::math::check_pos_definite("cholesky_decompose", "m", llt);
  stan::math::matrix_d m1_cpu = llt.matrixL();

  stan::math::matrix_cl<double> m1_pos_def_cl(m1_pos_def);
  stan::math::matrix_cl<double> m1_pos_def_cl_result
      = stan::math::cholesky_decompose(m1_pos_def_cl);
  stan::math::matrix_d m1_cl_cpu
      = stan::math::from_matrix_cl(m1_pos_def_cl_result);
  double max_error = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j <= i; j++) {
      double abs_err = std::fabs(m1_cpu(i, j) - m1_cl_cpu(i, j));
      double a = std::max(abs_err / m1_cpu(i, j), abs_err / m1_cl_cpu(i, j));
      max_error = std::max(max_error, a);
    }
  }
}
}  // namespace

TEST(MathMatrixOpenCLPrim, cholesky_decompose_small) {
  cholesky_decompose_test(10);
  cholesky_decompose_test(50);
  cholesky_decompose_test(100);
}

TEST(MathMatrixOpenCLPrim, cholesky_decompose_big) {
  cholesky_decompose_test(1251);
  cholesky_decompose_test(1704);
  cholesky_decompose_test(2000);
}
#endif
