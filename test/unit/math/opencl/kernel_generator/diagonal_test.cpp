#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using stan::math::diagonal;
using stan::math::matrix_cl;

TEST(KernelGenerator, diagonal_test) {
  MatrixXd m = MatrixXd::Random(3, 4);

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = diagonal(m_cl);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m.diagonal();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, diagonal_multiple_operations_test) {
  MatrixXd m = MatrixXd::Random(4, 3);

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = diagonal(m_cl * 2);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = (2 * m).diagonal();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, diagonal_multiple_operations_accept_lvalue_test) {
  MatrixXd m = MatrixXd::Random(3, 4);

  matrix_cl<double> m_cl(m);

  auto tmp = m_cl * 2;
  matrix_cl<double> res_cl = diagonal(tmp);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = (2 * m).diagonal();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, diagonal_lhs_test) {
  MatrixXd m = MatrixXd::Random(3, 4);
  VectorXd v = VectorXd::Random(3);

  matrix_cl<double> m_cl(m);
  matrix_cl<double> v_cl(v);

  diagonal(m_cl) = v_cl;
  MatrixXd res = stan::math::from_matrix_cl(m_cl);

  MatrixXd correct = m;
  correct.diagonal() = v;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, diagonal_of_a_block_test) {
  MatrixXd m = MatrixXd::Random(4, 4);

  matrix_cl<double> m_cl(m);

  diagonal(block_zero_based(m_cl, 0, 1, 3, 3))
      = diagonal(block_zero_based(m_cl, 1, 0, 3, 3));
  MatrixXd res = stan::math::from_matrix_cl(m_cl);

  MatrixXd correct = m;
  correct.diagonal(1) = correct.diagonal(-1);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, diagonal_lhs_add_assign_test) {
  MatrixXd m = MatrixXd::Random(3, 4);
  MatrixXd correct = m;
  correct.diagonal().array() += 1;

  matrix_cl<double> m_cl(m);

  diagonal(m_cl) += 1;
  MatrixXd res = stan::math::from_matrix_cl(m_cl);

  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
