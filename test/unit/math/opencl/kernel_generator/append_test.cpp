#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/prim/fun.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

TEST(KernelGenerator, append_row_errors) {
  using Eigen::MatrixXd;
  using stan::math::append_row;
  using stan::math::matrix_cl;

  MatrixXd m1(1, 1);
  m1 << 1;
  MatrixXd rv(1, 3);
  rv << 1, 2.5, 3;
  MatrixXd m2(2, 2);
  m2 << 1, 2.5, 3, 4;
  MatrixXd m3(3, 3);
  m3 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  matrix_cl<double> m3_cl(m3);
  matrix_cl<double> rv_cl(rv);

  EXPECT_NO_THROW(append_row(m3_cl, rv_cl));
  EXPECT_NO_THROW(append_row(m3_cl, stan::math::rowwise_broadcast(m1_cl)));
  EXPECT_NO_THROW(append_row(stan::math::rowwise_broadcast(m1_cl), m3_cl));
  EXPECT_THROW(append_row(m3_cl, stan::math::colwise_broadcast(rv_cl)),
               std::invalid_argument);
  EXPECT_THROW(append_row(stan::math::colwise_broadcast(rv_cl), m3_cl),
               std::invalid_argument);
  EXPECT_THROW(append_row(m2_cl, m3_cl), std::invalid_argument);
}

TEST(KernelGenerator, append_row_test) {
  using Eigen::MatrixXd;
  using stan::math::append_row;
  using stan::math::matrix_cl;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(1, 3);
  m2 << 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  matrix_cl<double> res_cl = append_row(m1_cl, m2_cl);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = append_row(m1, m2);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, append_row_multiple_operations_test) {
  using Eigen::MatrixXd;
  using stan::math::append_row;
  using stan::math::matrix_cl;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(1, 3);
  m2 << 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  matrix_cl<double> res_cl
      = append_row(append_row(m1_cl, m1_cl), append_row(m2_cl, m2_cl));
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = append_row(append_row(m1, m1), append_row(m2, m2));
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, append_row_multiple_operations_accept_lvalue_test) {
  using Eigen::MatrixXd;
  using stan::math::append_row;
  using stan::math::matrix_cl;
  MatrixXd m1(2, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(1, 3);
  m2 << 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto a = append_row(m1_cl, m1_cl);
  auto b = append_row(m2_cl, m2_cl);
  matrix_cl<double> res_cl = append_row(a, b);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = append_row(append_row(m1, m1), append_row(m2, m2));
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, append_col_errors) {
  using Eigen::MatrixXd;
  using stan::math::append_col;
  using stan::math::matrix_cl;

  MatrixXd m1(1, 1);
  m1 << 1;
  MatrixXd v(3, 1);
  v << 1, 2.5, 3;
  MatrixXd m2(2, 2);
  m2 << 1, 2.5, 3, 4;
  MatrixXd m3(3, 3);
  m3 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  matrix_cl<double> m3_cl(m3);
  matrix_cl<double> v_cl(v);

  EXPECT_NO_THROW(append_col(m3_cl, v_cl));
  EXPECT_NO_THROW(append_col(m3_cl, stan::math::colwise_broadcast(m1_cl)));
  EXPECT_NO_THROW(append_col(stan::math::colwise_broadcast(m1_cl), m3_cl));
  EXPECT_THROW(append_col(m3_cl, stan::math::rowwise_broadcast(v_cl)),
               std::invalid_argument);
  EXPECT_THROW(append_col(stan::math::rowwise_broadcast(v_cl), m3_cl),
               std::invalid_argument);
  EXPECT_THROW(append_col(m2_cl, m3_cl), std::invalid_argument);
}

TEST(KernelGenerator, append_col_test) {
  using Eigen::MatrixXd;
  using stan::math::append_col;
  using stan::math::matrix_cl;
  MatrixXd m1(3, 2);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(3, 1);
  m2 << 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  matrix_cl<double> res_cl = append_col(m1_cl, m2_cl);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = append_col(m1, m2);
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, append_col_multiple_operations_test) {
  using Eigen::MatrixXd;
  using stan::math::append_col;
  using stan::math::matrix_cl;
  MatrixXd m1(3, 2);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(3, 1);
  m2 << 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  matrix_cl<double> res_cl
      = append_col(append_col(m1_cl, m1_cl), append_col(m2_cl, m2_cl));
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = append_col(append_col(m1, m1), append_col(m2, m2));
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, append_col_multiple_operations_accept_lvalue_test) {
  using Eigen::MatrixXd;
  using stan::math::append_col;
  using stan::math::matrix_cl;
  MatrixXd m1(3, 2);
  m1 << 1, 2.5, 3, 4, 5, 6.3;
  MatrixXd m2(3, 1);
  m2 << 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto a = append_col(m1_cl, m1_cl);
  auto b = append_col(m2_cl, m2_cl);
  matrix_cl<double> res_cl = append_col(a, b);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = append_col(append_col(m1, m1), append_col(m2, m2));
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#endif
