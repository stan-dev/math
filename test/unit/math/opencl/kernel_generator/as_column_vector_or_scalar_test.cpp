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

TEST(KernelGenerator, as_column_vector_or_scalar_errors) {
  using stan::math::as_column_vector_or_scalar;
  stan::math::matrix_cl<double> m(7, 1);
  stan::math::matrix_cl<double> n(1, 7);
  stan::math::matrix_cl<double> j(7, 7);

  auto block_m = stan::math::block_zero_based(m, 0, 0, 7, 1);
  auto block_n = stan::math::block_zero_based(n, 0, 0, 1, 7);

  EXPECT_THROW(as_column_vector_or_scalar(j), std::invalid_argument);
  EXPECT_THROW(block_n = as_column_vector_or_scalar(n), std::invalid_argument);
  EXPECT_THROW(block_n = as_column_vector_or_scalar(m), std::invalid_argument);
  EXPECT_NO_THROW(block_m = as_column_vector_or_scalar(m));
  EXPECT_NO_THROW(as_column_vector_or_scalar(block_m) = m);
  EXPECT_NO_THROW(as_column_vector_or_scalar(block_m)
                  = as_column_vector_or_scalar(m));
  EXPECT_NO_THROW(block_m = as_column_vector_or_scalar(n));
  EXPECT_NO_THROW(as_column_vector_or_scalar(block_n) = m);
  EXPECT_NO_THROW(as_column_vector_or_scalar(block_n)
                  = as_column_vector_or_scalar(n));
  auto a = as_column_vector_or_scalar(m);
  EXPECT_NO_THROW(a = a);
  EXPECT_NO_THROW(a = a + 1);
}

TEST(KernelGenerator, as_column_vector_or_scalar_vector_test) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::matrix_cl;
  VectorXd m(6, 1);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp = stan::math::as_column_vector_or_scalar(m_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::as_column_vector_or_scalar(m);
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, as_column_vector_or_scalar_row_vector_test) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using stan::math::matrix_cl;
  RowVectorXd m(6);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp = stan::math::as_column_vector_or_scalar(m_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::as_column_vector_or_scalar(m);
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, double_as_column_vector_or_scalar_test) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::matrix_cl;
  VectorXd m(6);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp = stan::math::as_column_vector_or_scalar(
      stan::math::as_column_vector_or_scalar(m_cl));
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m;
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, double_as_column_vector_or_scalar_accepts_lvalue_test) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using stan::math::matrix_cl;
  RowVectorXd m(6);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);
  auto tmp2 = stan::math::as_column_vector_or_scalar(m_cl);
  auto tmp = stan::math::as_column_vector_or_scalar(tmp2);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = stan::math::as_column_vector_or_scalar(m);
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, as_column_vector_or_scalar_vector_block_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(5, 5);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      21, 22, 23, 24, 25;

  matrix_cl<double> m_cl(m);
  auto tmp2 = stan::math::block_zero_based(m_cl, 1, 2, 3, 1);
  auto tmp = stan::math::as_column_vector_or_scalar(tmp2);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct
      = stan::math::as_column_vector_or_scalar(m.col(2).segment(1, 3));
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, as_column_vector_or_scalar_row_vector_block_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(5, 5);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
      21, 22, 23, 24, 25;

  matrix_cl<double> m_cl(m);
  auto tmp2 = stan::math::block_zero_based(m_cl, 2, 1, 1, 3);
  auto tmp = stan::math::as_column_vector_or_scalar(tmp2);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct
      = stan::math::as_column_vector_or_scalar(m.row(2).segment(1, 3));
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(KernelGenerator, lhs_as_column_vector_or_scalar_vector_plus_eq_test) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::as_column_vector_or_scalar;
  using stan::math::matrix_cl;
  VectorXd m1(6);
  m1 << 1, 2, 3, 4, 5, 6;
  VectorXd correct = m1.array() + 1;

  matrix_cl<double> m1_cl(m1);

  as_column_vector_or_scalar(m1_cl) += 1;

  MatrixXd res = stan::math::from_matrix_cl(m1_cl);

  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, lhs_as_column_vector_or_scalar_row_vector_plus_eq_test) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using stan::math::as_column_vector_or_scalar;
  using stan::math::matrix_cl;
  RowVectorXd m1(6);
  m1 << 1, 2, 3, 4, 5, 6;
  RowVectorXd correct = m1.array() + 1;

  matrix_cl<double> m1_cl(m1);

  as_column_vector_or_scalar(m1_cl) += 1;

  MatrixXd res = stan::math::from_matrix_cl(m1_cl);

  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
