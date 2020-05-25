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

TEST(KernelGenerator, block_errors) {
  using stan::math::block;

  matrix_cl<double> m(7, 9);

  EXPECT_NO_THROW(block(m, 0, 0, 1, 1));
  EXPECT_NO_THROW(block(m, 6, 8, 1, 1));
  EXPECT_NO_THROW(block(m, 0, 0, 7, 9));
  EXPECT_THROW(block(m, 0, 0, 8, 1), std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, 1, 10), std::invalid_argument);
  EXPECT_THROW(block(m, 6, -1, 1, 1), std::invalid_argument);
  EXPECT_THROW(block(m, -1, 5, 1, 1), std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, -1, 1), std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, 1, -1), std::invalid_argument);
  EXPECT_THROW(block(m, 0, 9, 0, 1), std::invalid_argument);

  EXPECT_NO_THROW(block(m, 0, 0, 7, 9) = m);
  EXPECT_THROW(block(m, 0, 0, 7, 8) = m, std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, 6, 9) = m, std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, 7, 8) = m, std::invalid_argument);
  EXPECT_THROW(block(m, 0, 0, 6, 9) = m, std::invalid_argument);
}

TEST(KernelGenerator, block_test) {
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

TEST(KernelGenerator, block_multiple_operations_test) {
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

TEST(KernelGenerator, block_multiple_operations_accept_lvalue_test) {
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

TEST(KernelGenerator, lhs_block_test) {
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

TEST(KernelGenerator, block_to_lhs_block_test) {
  using stan::math::block;
  MatrixXd m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  MatrixXd m2 = MatrixXd::Constant(5, 7, 9);

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  block(m2_cl, 1, 1, 2, 3) = block(m1_cl, 0, 0, 2, 3);

  MatrixXd res = stan::math::from_matrix_cl(m2_cl);

  MatrixXd correct = m2;
  correct.block(1, 1, 2, 3) = m1;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, two_blocks_of_same_expression) {
  using stan::math::block;
  MatrixXd m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  matrix_cl<double> m_cl(m);

  auto tmp = m_cl + 1;
  auto tmp2 = block(tmp, 0, 0, 2, 2) + block(tmp, 0, 1, 2, 2);

  matrix_cl<double> res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = (m.block(0, 0, 2, 2) + m.block(0, 1, 2, 2)).array() + 2;

  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, block_view_test) {
  using stan::math::block;
  matrix_cl<double> m(4, 4, stan::math::matrix_cl_view::Diagonal);
  matrix_cl<double> res = block(m, 0, 0, 2, 2);
  EXPECT_EQ(res.view(), stan::math::matrix_cl_view::Diagonal);
  res = block(m, 1, 0, 2, 2);
  EXPECT_EQ(res.view(), stan::math::matrix_cl_view::Upper);
  res = block(m, 0, 1, 2, 2);
  EXPECT_EQ(res.view(), stan::math::matrix_cl_view::Lower);
  res = block(m, 0, 2, 2, 2);
  EXPECT_EQ(res.view(), stan::math::matrix_cl_view::Diagonal);
  res = block(m, 2, 0, 2, 2);
  EXPECT_EQ(res.view(), stan::math::matrix_cl_view::Diagonal);

  res = block(cos(m), 1, 0, 2, 2);
  EXPECT_EQ(res.view(), stan::math::matrix_cl_view::Entire);
}

#endif
