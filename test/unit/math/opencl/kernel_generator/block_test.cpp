#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  EXPECT_EQ(A.rows(), B.rows());        \
  EXPECT_EQ(A.cols(), B.cols());        \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixCL, block_errors) {
  using stan::math::block;

  matrix_cl<double> m(7, 9);

  EXPECT_NO_THROW(block(m, 0, 0, 1, 1));
  EXPECT_NO_THROW(block(m, 6, 8, 1, 1));
  EXPECT_NO_THROW(block(m, 0, 0, 7, 9));
  EXPECT_THROW(block(m, 0, 0, 8, 1), std::domain_error);
  EXPECT_THROW(block(m, 0, 0, 1, 10), std::domain_error);
  EXPECT_THROW(block(m, 8, 0, 1, 1), std::domain_error);
  EXPECT_THROW(block(m, 0, 9, 1, 1), std::domain_error);

  EXPECT_NO_THROW(block(m, 0, 0, 7, 9) = m);
  EXPECT_THROW(block(m, 0, 0, 7, 8) = m, std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, 6, 9) = m, std::invalid_argument);
}

TEST(MathMatrixCL, block_test) {
  using stan::math::block;
  std::string kernel_filename = "block.cl";
  MatrixXd m = MatrixXd::Random(7, 9);

  matrix_cl<double> m_cl(m);

  auto tmp = block(m_cl, 2, 4, 3, 5);
  matrix_cl<double> res_cl;
  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m.block(2, 4, 3, 5);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, block_multiple_operations_test) {
  using stan::math::block;
  MatrixXd m = MatrixXd::Random(7, 9);

  matrix_cl<double> m_cl(m);

  auto tmp = block(block(m_cl, 1, 1, 5, 5), 2, 0, 2, 3);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd tmp_block = m.block(1, 1, 5, 5);
  MatrixXd correct = tmp_block.block(2, 0, 2, 3);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, block_multiple_operations_accept_lvalue_test) {
  using stan::math::block;
  MatrixXd m = MatrixXd::Random(7, 9);

  matrix_cl<double> m_cl(m);

  auto tmp2 = block(m_cl, 1, 1, 5, 5);
  auto tmp = block(tmp2, 2, 0, 2, 3);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd tmp_block = m.block(1, 1, 5, 5);
  MatrixXd correct = tmp_block.block(2, 0, 2, 3);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, lhs_block_test) {
  using stan::math::block;
  MatrixXd m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  MatrixXd m2 = MatrixXd::Constant(5, 7, 9);

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  block(m2_cl, 1, 1, 2, 3) = m1_cl;

  MatrixXd res = stan::math::from_matrix_cl(m2_cl);

  MatrixXd correct = m2;
  correct.block(1, 1, 2, 3) = m1;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
