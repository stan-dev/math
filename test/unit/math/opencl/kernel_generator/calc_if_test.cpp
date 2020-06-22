#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, is_without_output_test) {
  EXPECT_TRUE((stan::math::is_without_output<
               stan::math::calc_if_<false, stan::math::scalar_<int>>>::value));
  EXPECT_TRUE((stan::math::is_without_output<const stan::math::calc_if_<
                   false, stan::math::scalar_<int>>&>::value));
  EXPECT_FALSE((stan::math::is_without_output<
                stan::math::calc_if_<true, stan::math::scalar_<int>>>::value));
  EXPECT_FALSE(
      (stan::math::is_without_output<stan::math::scalar_<int>>::value));
}

TEST(KernelGenerator, calc_if_pass_test) {
  double eps = 0;
  MatrixXd m(2, 2);
  m << 1, 2, 3, 4;
  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl = stan::math::calc_if<true>(m_cl);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  EXPECT_MATRIX_NEAR(res, m, eps);
}

TEST(KernelGenerator, calc_if_multi_result_pass_test) {
  double eps = 0;
  MatrixXd m(2, 2);
  m << 1, 2, 3, 4;
  matrix_cl<double> m_cl(m);

  auto tmp = stan::math::calc_if<true>(m_cl);
  matrix_cl<double> res_cl;
  stan::math::results(res_cl) = stan::math::expressions(tmp);
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  EXPECT_MATRIX_NEAR(res, m, eps);
}

TEST(KernelGenerator, calc_if_multi_result_no_output_test) {
  double eps = 0;
  MatrixXd m(2, 2);
  m << 1, 2, 3, 4;
  matrix_cl<double> m_cl(m);

  auto tmp = stan::math::calc_if<false>(m_cl);
  matrix_cl<double> res_cl;
  stan::math::results(res_cl) = stan::math::expressions(tmp);
  EXPECT_EQ(res_cl.size(), 0);
}

TEST(KernelGenerator, calc_if_multi_result_one_output_test) {
  double eps = 0;
  MatrixXd m(2, 2);
  m << 1, 2, 3, 4;
  matrix_cl<double> m_cl(m);

  matrix_cl<double> res1_cl;
  matrix_cl<double> res2_cl;
  stan::math::results(res1_cl, res2_cl) = stan::math::expressions(
      stan::math::calc_if<false>(m_cl + 1), stan::math::calc_if<true>(m_cl));
  EXPECT_EQ(res1_cl.size(), 0);

  MatrixXd res = stan::math::from_matrix_cl(res2_cl);
  EXPECT_MATRIX_NEAR(res, m, eps);
}

TEST(KernelGenerator, calc_if_colwise_min_test) {
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl
      = stan::math::calc_if<true>(stan::math::colwise_min(m_cl));
  MatrixXd raw_res = stan::math::from_matrix_cl(res_cl);
  EXPECT_GE(m.rows(), raw_res.rows());
  MatrixXd res = raw_res.colwise().minCoeff();
  MatrixXd correct = m.colwise().minCoeff();
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

#endif
