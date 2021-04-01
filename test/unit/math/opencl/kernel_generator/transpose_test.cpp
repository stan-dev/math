#ifdef STAN_OPENCL
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(KernelGenerator, transpose_errors) {
  using stan::math::transpose;
  stan::math::matrix_cl<double> m(7, 9);
  stan::math::matrix_cl<double> n(7, 9);

  EXPECT_THROW(block(n, 0, 0, 7, 9) = transpose(m), std::invalid_argument);
  EXPECT_THROW(transpose(block(n, 0, 0, 7, 9)) = m, std::invalid_argument);
  EXPECT_NO_THROW(n = transpose(m));
  EXPECT_NO_THROW(transpose(m) = n);
  EXPECT_NO_THROW(transpose(m) = transpose(m));
  auto a = transpose(m);
  EXPECT_NO_THROW(a = a);
  EXPECT_NO_THROW(a = a + 1);
}

TEST(KernelGenerator, transpose_rvalue_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  auto tmp = stan::math::transpose(stan::math::to_matrix_cl(m));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::transpose(m);
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, transpose_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp = stan::math::transpose(m_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::transpose(m);
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, transpose_triangular_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m, stan::math::matrix_cl_view::Upper);
  auto tmp = stan::math::transpose(m_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::transpose(m).triangularView<Eigen::Lower>();
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Lower, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, double_transpose_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp = stan::math::transpose(stan::math::transpose(m_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m;
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, double_transpose_triangular_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m, stan::math::matrix_cl_view::Upper);
  auto tmp = stan::math::transpose(stan::math::transpose(m_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m.triangularView<Eigen::Upper>();
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Upper, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, double_transpose_accepts_lvalue_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp2 = stan::math::transpose(m_cl);
  auto tmp = stan::math::transpose(tmp2);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m;
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, transpose_block_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(5, 5);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      21, 22, 23, 24, 25;

  matrix_cl<double> m_cl(m);
  auto tmp2 = stan::math::block_zero_based(m_cl, 2, 1, 2, 3);
  auto tmp = stan::math::transpose(tmp2);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m.block(2, 1, 2, 3).transpose();
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, block_of_transpose_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(5, 5);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      21, 22, 23, 24, 25;

  matrix_cl<double> m_cl(m);
  auto tmp = stan::math::transpose(m_cl);
  auto tmp2 = stan::math::block_zero_based(tmp, 2, 1, 2, 3);
  matrix_cl<double> res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m.transpose().block(2, 1, 2, 3);
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, a_plus_a_transpose_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  std::string kernel_filename = "a+aT.cl";
  MatrixXd m(3, 3);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9;

  matrix_cl<double> m_cl(m);
  auto tmp2 = m_cl + 1;
  auto tmp = tmp2 + stan::math::transpose(tmp2);
  matrix_cl<double> res_cl;

  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m.array() + 2 + m.array().transpose();
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_EQ(stan::math::matrix_cl_view::Entire, res_cl.view());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, lhs_transpose_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  using stan::math::transpose;
  MatrixXd m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  MatrixXd m2 = MatrixXd::Constant(5, 3, 2);
  MatrixXd correct = m2;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  transpose(m2_cl) = m1_cl;

  MatrixXd res = stan::math::from_matrix_cl(m2_cl);

  correct.transpose() = m1;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, transpose_to_lhs_transpose_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  using stan::math::transpose;
  MatrixXd m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  MatrixXd m2 = MatrixXd::Constant(5, 3, 2);
  MatrixXd correct = m2;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  transpose(m2_cl) = transpose(m1_cl);

  MatrixXd res = stan::math::from_matrix_cl(m2_cl);

  correct.transpose() = m1.transpose();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, lhs_transpose_plus_eq_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  using stan::math::transpose;
  MatrixXd m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  MatrixXd correct = m1.array() + 1;

  matrix_cl<double> m1_cl(m1);

  transpose(m1_cl) += 1;

  MatrixXd res = stan::math::from_matrix_cl(m1_cl);

  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
