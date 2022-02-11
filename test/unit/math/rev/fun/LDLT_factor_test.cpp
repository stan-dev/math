#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(AgradRevMatrix, LDLT_factor_default_constructor) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;
  stan::math::var_value<Eigen::MatrixXd> A_vm = A.val();

  EXPECT_NO_THROW(stan::math::make_ldlt_factor(Eigen::MatrixXd()));
  EXPECT_NO_THROW(stan::math::make_ldlt_factor(A));
  EXPECT_NO_THROW(stan::math::make_ldlt_factor(A_vm));

  auto ldlt_A = stan::math::make_ldlt_factor(A);
  EXPECT_EQ(2U, ldlt_A.matrix().rows());
  EXPECT_EQ(2U, ldlt_A.matrix().cols());

  auto ldlt_A_vm = stan::math::make_ldlt_factor(A_vm);
  EXPECT_EQ(2U, ldlt_A_vm.matrix().rows());
  EXPECT_EQ(2U, ldlt_A_vm.matrix().cols());

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, solve) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;
  stan::math::var_value<Eigen::MatrixXd> A_vm = A.val();
  Eigen::Matrix<double, -1, -1> B(2, 2);
  B(0, 0) = 3;
  B(0, 1) = 1;
  B(1, 0) = 2;
  B(1, 1) = 2;

  Eigen::Matrix<double, -1, -1> expected_solve(2, 2);
  expected_solve(0, 0) = 1.333333333;
  expected_solve(0, 1) = 0.0;
  expected_solve(1, 0) = 0.333333333;
  expected_solve(1, 1) = 1.0;
  Eigen::Matrix<double, -1, -1> solve;

  auto ldlt_A = stan::math::make_ldlt_factor(A);
  EXPECT_NO_THROW(solve = ldlt_A.ldlt().solve(B));

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      EXPECT_FLOAT_EQ(expected_solve(i, j), solve(i, j));

  auto ldlt_A_vm = stan::math::make_ldlt_factor(A_vm);
  EXPECT_NO_THROW(solve = ldlt_A_vm.ldlt().solve(B));

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      EXPECT_FLOAT_EQ(expected_solve(i, j), solve(i, j));

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, matrix) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;
  auto ldlt_A = stan::math::make_ldlt_factor(A);

  stan::math::var_value<Eigen::MatrixXd> A_vm = A.val();
  auto ldlt_A_vm = stan::math::make_ldlt_factor(A_vm);

  Eigen::Matrix<double, -1, -1> A_double(2, 2);
  A_double << 5, 1, 1, 5;
  auto ldlt_double = stan::math::make_ldlt_factor(A_double);

  EXPECT_MATRIX_EQ(A_double, ldlt_double.matrix());
  EXPECT_MATRIX_EQ(A.val(), ldlt_A.matrix().val());
  EXPECT_MATRIX_EQ(A_vm.val(), ldlt_A_vm.matrix().val());

  stan::math::recover_memory();
}
