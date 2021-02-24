#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>

TEST(KernelGenerator, indexing_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;

  std::string kernel_filename = "indexing.cl";
  MatrixXd m = MatrixXd::Random(7, 9);

  Eigen::MatrixXi col_idx(3, 2);
  col_idx << 2, 5, 0, 1, 6, 3;
  Eigen::MatrixXi row_idx(3, 2);
  row_idx << 1, 2, 3, 4, 2, 0;

  matrix_cl<double> m_cl(m);
  matrix_cl<int> row_idx_cl(row_idx);
  matrix_cl<int> col_idx_cl(col_idx);

  auto tmp = stan::math::indexing(m_cl, row_idx_cl, col_idx_cl);

  matrix_cl<double> res_cl;
  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct(3, 2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      correct(i, j) = m(row_idx(i, j), col_idx(i, j));
    }
  }
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(KernelGenerator, indexing_multiple_operations_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m = MatrixXd::Random(7, 9);

  Eigen::MatrixXi col_idx(3, 2);
  col_idx << 2, 5, 0, 1, 6, 3;
  Eigen::MatrixXi row_idx(3, 2);
  row_idx << 1, 2, 3, 4, 2, 0;

  Eigen::MatrixXi col_idx2(2, 2);
  col_idx2 << 0, 0, 1, 1;
  Eigen::MatrixXi row_idx2(2, 2);
  row_idx2 << 0, 1, 1, 2;

  matrix_cl<double> m_cl(m);
  matrix_cl<int> row_idx_cl(row_idx);
  matrix_cl<int> col_idx_cl(col_idx);
  matrix_cl<int> row_idx2_cl(row_idx2);
  matrix_cl<int> col_idx2_cl(col_idx2);

  matrix_cl<double> res_cl = stan::math::indexing(
      stan::math::indexing(m_cl, row_idx_cl, col_idx_cl),
      stan::math::indexing(row_idx2_cl, col_idx2_cl, col_idx2_cl), col_idx2_cl);

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd tmp(3, 2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      tmp(i, j) = m(row_idx(i, j), col_idx(i, j));
    }
  }

  MatrixXd correct(2, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      int a = col_idx2(i, j);
      correct(i, j) = tmp(row_idx2(a, a), a);
    }
  }
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(KernelGenerator, indexing_lhs_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m = MatrixXd::Zero(7, 9);
  MatrixXd m2 = MatrixXd::Random(3, 2);

  Eigen::MatrixXi col_idx(3, 2);
  col_idx << 2, 5, 0, 1, 6, 3;
  Eigen::MatrixXi row_idx(3, 2);
  row_idx << 1, 2, 3, 4, 2, 0;

  matrix_cl<double> m_cl(m);
  matrix_cl<double> m2_cl(m2);
  matrix_cl<int> row_idx_cl(row_idx);
  matrix_cl<int> col_idx_cl(col_idx);

  stan::math::indexing(m_cl, row_idx_cl, col_idx_cl) = m2_cl;

  MatrixXd res = stan::math::from_matrix_cl(m_cl);

  MatrixXd correct = m;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      correct(row_idx(i, j), col_idx(i, j)) = m2(i, j);
    }
  }
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(KernelGenerator, indexing_repeat_lhs_rhs_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m = MatrixXd::Zero(7, 9);
  MatrixXd correct = m;

  Eigen::MatrixXi col_idx(3, 2);
  col_idx << 2, 5, 0, 1, 6, 3;
  Eigen::MatrixXi row_idx(3, 2);
  row_idx << 1, 2, 3, 4, 2, 0;

  matrix_cl<double> m_cl(m);
  matrix_cl<int> row_idx_cl(row_idx);
  matrix_cl<int> col_idx_cl(col_idx);

  auto b = stan::math::indexing(m_cl, row_idx_cl, col_idx_cl);

  b = b + 1;
  MatrixXd res = stan::math::from_matrix_cl(m_cl);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      correct(row_idx(i, j), col_idx(i, j)) += 1;
    }
  }
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(KernelGenerator, indexing_lhs_plus_assign_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m = MatrixXd::Zero(7, 9);
  MatrixXd correct = m;

  Eigen::MatrixXi col_idx(3, 2);
  col_idx << 2, 5, 0, 1, 6, 3;
  Eigen::MatrixXi row_idx(3, 2);
  row_idx << 1, 2, 3, 4, 2, 0;

  matrix_cl<double> m_cl(m);
  matrix_cl<int> row_idx_cl(row_idx);
  matrix_cl<int> col_idx_cl(col_idx);

  auto b = stan::math::indexing(m_cl, row_idx_cl, col_idx_cl);

  b += 1;
  MatrixXd res = stan::math::from_matrix_cl(m_cl);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      correct(row_idx(i, j), col_idx(i, j)) += 1;
    }
  }
  EXPECT_MATRIX_EQ(res, correct);
}

#endif
