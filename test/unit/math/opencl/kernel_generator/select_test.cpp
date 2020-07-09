#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, select_errors) {
  using stan::math::select;

  matrix_cl<double> m_double(2, 3);
  matrix_cl<double> m_rows(3, 3);
  matrix_cl<double> m_cols(2, 4);
  matrix_cl<int> m_int(2, 3);
  matrix_cl<bool> m_bool(2, 3);

  EXPECT_NO_THROW(select(m_bool, m_double, m_double));
  EXPECT_NO_THROW(select(m_bool, m_double, m_int));
  EXPECT_NO_THROW(select(m_bool, m_int, m_double));
  EXPECT_NO_THROW(select(m_int, m_double, m_double));
  EXPECT_THROW(select(m_rows, m_double, m_double), std::invalid_argument);
  EXPECT_THROW(select(m_cols, m_double, m_double), std::invalid_argument);
  EXPECT_THROW(select(m_bool, m_rows, m_double), std::invalid_argument);
  EXPECT_THROW(select(m_bool, m_cols, m_double), std::invalid_argument);
  EXPECT_THROW(select(m_bool, m_double, m_rows), std::invalid_argument);
  EXPECT_THROW(select(m_bool, m_double, m_cols), std::invalid_argument);
}

TEST(KernelGenerator, select_zero_size) {
  matrix_cl<double> m_bool_zero(0, 0);
  matrix_cl<double> m_double_zero(0, 0);
  matrix_cl<int> m_int_zero(0, 0);
  matrix_cl<int> m_int_zero_rows(0, 1);

  EXPECT_NO_THROW(select(m_bool_zero, m_double_zero, m_double_zero));
  EXPECT_NO_THROW(select(m_bool_zero, m_double_zero, m_int_zero));
  EXPECT_NO_THROW(select(m_bool_zero, m_int_zero, m_double_zero));
  EXPECT_THROW(select(m_bool_zero, m_double_zero, m_int_zero_rows),
               std::invalid_argument);
}

TEST(KernelGenerator, select_test) {
  std::string kernel_filename = "select.cl";
  using stan::math::select;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(2, 3);
  m2 << 10, 100, 1000, 30.5, 12.1, 0.3;
  Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_cond(2, 3);
  m_cond << true, false, false, true, true, false;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  matrix_cl<bool> m_cond_cl(m_cond);

  auto tmp = select(m_cond_cl, m1_cl, m2_cl);

  matrix_cl<double> res_cl;
  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m_cond.select(m1, m2);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, multiple_operations) {
  using stan::math::select;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXi m2(2, 3);
  m2 << 10, 100, 1000, 30, 12, 0;
  MatrixXd m3(2, 3);
  m3 << 1.4, 12.5, 6.3, 34.4, -5, -6.3;
  Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_cond(2, 3);
  m_cond << true, false, false, true, true, false;
  Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_cond2(2, 3);
  m_cond2 << false, false, true, true, false, true;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);
  matrix_cl<double> m3_cl(m3);
  matrix_cl<bool> m_cond_cl(m_cond);
  matrix_cl<bool> m_cond2_cl(m_cond2);
  auto tmp = select(m_cond_cl, m1_cl, select(m_cond2_cl, m2_cl, m3_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m_cond.select(m1, m_cond2.select(m2.cast<double>(), m3));
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, multiple_operations_size_zero) {
  using stan::math::select;

  matrix_cl<double> m1_cl(0, 0);
  matrix_cl<int> m2_cl(0, 0);
  matrix_cl<double> m3_cl(0, 0);
  matrix_cl<bool> m_cond_cl(0, 0);
  matrix_cl<bool> m_cond2_cl(0, 0);
  auto tmp = select(m_cond_cl, m1_cl, select(m_cond2_cl, m2_cl, m3_cl));
  matrix_cl<double> res_cl = tmp;

  EXPECT_EQ(0, res_cl.rows());
  EXPECT_EQ(0, res_cl.cols());
}

TEST(KernelGenerator, multiple_operations_accepts_lvalue) {
  using stan::math::select;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXi m2(2, 3);
  m2 << 10, 100, 1000, 30, 12, 0;
  MatrixXd m3(2, 3);
  m3 << 1.4, 12.5, 6.3, 34.4, -5, -6.3;
  Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_cond(2, 3);
  m_cond << true, false, false, true, true, false;
  Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_cond2(2, 3);
  m_cond2 << false, false, true, true, false, true;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);
  matrix_cl<double> m3_cl(m3);
  matrix_cl<bool> m_cond_cl(m_cond);
  matrix_cl<bool> m_cond2_cl(m_cond2);
  auto tmp2 = select(m_cond2_cl, m2_cl, m3_cl);
  auto tmp = select(m_cond_cl, m1_cl, tmp2);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m_cond.select(m1, m_cond2.select(m2.cast<double>(), m3));
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
