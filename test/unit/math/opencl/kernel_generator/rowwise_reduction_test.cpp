#ifdef STAN_OPENCL

#include <stan/math.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::MatrixXd;
using stan::math::matrix_cl;

TEST(KernelGenerator, rowwise_sum_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().sum();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_prod_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_prod(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().prod();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_min_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_min(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().minCoeff();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_max_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_max(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().maxCoeff();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_sum_triangular_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m, stan::math::matrix_cl_view::Lower);

  MatrixXd m12 = m;
  m12.triangularView<Eigen::StrictlyUpper>() = MatrixXd::Constant(3, 2, 0);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m12.rowwise().sum();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);

  m_cl.view(stan::math::matrix_cl_view::Upper);

  MatrixXd m34 = m;
  m34.triangularView<Eigen::StrictlyLower>() = MatrixXd::Constant(3, 2, 0);

  matrix_cl<double> res3_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res3 = stan::math::from_matrix_cl(res3_cl);
  MatrixXd correct3 = m34.rowwise().sum();
  EXPECT_EQ(correct3.rows(), res3.rows());
  EXPECT_EQ(correct3.cols(), res3.cols());
  EXPECT_MATRIX_NEAR(correct3, res3, 1e-9);
}

TEST(KernelGenerator, rowwise_min_triangular_test) {
  MatrixXd m(3, 2);
  m << -1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m, stan::math::matrix_cl_view::Lower);

  MatrixXd m12 = m;
  m12.triangularView<Eigen::StrictlyUpper>() = MatrixXd::Constant(3, 2, 0);

  matrix_cl<double> res1_cl = stan::math::rowwise_min(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m12.rowwise().minCoeff();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);

  m_cl.view(stan::math::matrix_cl_view::Upper);

  MatrixXd m34 = m;
  m34.triangularView<Eigen::StrictlyLower>() = MatrixXd::Constant(3, 2, 0);

  matrix_cl<double> res3_cl = stan::math::rowwise_min(m_cl);
  MatrixXd res3 = stan::math::from_matrix_cl(res3_cl);
  MatrixXd correct3 = m34.rowwise().minCoeff();
  EXPECT_EQ(correct3.rows(), res3.rows());
  EXPECT_EQ(correct3.cols(), res3.cols());
  EXPECT_MATRIX_NEAR(correct3, res3, 1e-9);
}

TEST(KernelGenerator, rowwise_sum_one_col_test) {
  MatrixXd m(3, 1);
  m << 1.1, 1.2, 1.3;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().sum();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_sum_one_row_test) {
  MatrixXd m(1, 3);
  m << 1.1, 1.2, 1.3;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().sum();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_sum_zero_rows_test) {
  matrix_cl<double> m_cl(0, 3);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  EXPECT_EQ(0, res1_cl.rows());
  EXPECT_EQ(1, res1_cl.cols());
}

TEST(KernelGenerator, rowwise_sum_zero_cols_test) {
  MatrixXd m(3, 0);

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().sum();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

TEST(KernelGenerator, rowwise_sum_zero_rows_and_cols_test) {
  matrix_cl<double> m_cl(0, 0);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  EXPECT_EQ(0, res1_cl.rows());
  EXPECT_EQ(1, res1_cl.cols());
}

TEST(KernelGenerator, rowwise_sum_big_test) {
  MatrixXd m = MatrixXd::Random(235, 135);

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl = stan::math::rowwise_sum(m_cl);
  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd correct1 = m.rowwise().sum();
  EXPECT_EQ(correct1.rows(), res1.rows());
  EXPECT_EQ(correct1.cols(), res1.cols());
  EXPECT_MATRIX_NEAR(correct1, res1, 1e-9);
}

#endif
