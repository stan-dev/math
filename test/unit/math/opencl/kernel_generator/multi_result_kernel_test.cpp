#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <string>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, multi_result_kernel_errors) {
  matrix_cl<double> m1_cl(3, 3);
  matrix_cl<double> m2_cl(3, 3);
  matrix_cl<double> m3_cl(4, 3);
  matrix_cl<double> m4_cl(4, 3);

  matrix_cl<double> res1_cl;
  matrix_cl<double> res2_cl;
  EXPECT_NO_THROW(stan::math::results() = stan::math::expressions());
  // mismatch in size between different (expression, result) pairs
  EXPECT_THROW(stan::math::results(res1_cl, res2_cl)
               = stan::math::expressions(m1_cl - m2_cl, m3_cl - m4_cl),
               std::invalid_argument);
  // mismatch in size between an expression and a result that can not be resized
  auto block1 = stan::math::block_zero_based(m1_cl, 0, 0, 3, 3);
  EXPECT_THROW(stan::math::results(res1_cl, block1)
               = stan::math::expressions(m3_cl + m4_cl, m3_cl - m4_cl),
               std::invalid_argument);
}

TEST(KernelGenerator, multi_result_kernel) {
  std::string kernel_filename = "sum_dif.cl";
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  matrix_cl<double> sum_cl;
  matrix_cl<double> diff_cl(3, 3);

  auto sum = m1_cl + m2_cl;
  auto res = stan::math::results(sum_cl, diff_cl);
  auto exprs = stan::math::expressions(sum, m1_cl - m2_cl);

  std::string kernel_src = res.get_kernel_source_for_evaluating(exprs);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res = exprs;

  MatrixXd res_sum = stan::math::from_matrix_cl(sum_cl);
  MatrixXd res_diff = stan::math::from_matrix_cl(diff_cl);

  MatrixXd correct_sum = m1 + m2;
  MatrixXd correct_diff = m1 - m2;

  EXPECT_MATRIX_NEAR(res_sum, correct_sum, 1e-9);
  EXPECT_MATRIX_NEAR(res_diff, correct_diff, 1e-9);
}

TEST(KernelGenerator, multi_result_kernel_add_assign) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;
  MatrixXd m3(3, 3);
  m3 << 10, 3, 1000, 2, 6, -142, 7, 6, 8;
  MatrixXd m4(3, 3);
  m4 << 6, 7, 8, 9, 1, 2, 3, 4, 5;

  MatrixXd m5 = m1, m6 = m2, m7 = m3, m8 = m4;

  matrix_cl<double> m1_cl(m5);
  matrix_cl<double> m2_cl(m6);
  matrix_cl<double> sum_cl(m7);
  matrix_cl<double> diff_cl(m8);

  auto sum = m1_cl + m2_cl;
  auto res = stan::math::results(sum_cl, diff_cl);
  auto exprs = stan::math::expressions(sum, m1_cl - m2_cl);

  res += exprs;

  MatrixXd res_sum = stan::math::from_matrix_cl(sum_cl);
  MatrixXd res_diff = stan::math::from_matrix_cl(diff_cl);

  m3 += m1 + m2;
  m4 += m1 - m2;

  EXPECT_MATRIX_NEAR(res_sum, m3, 1e-9);
  EXPECT_MATRIX_NEAR(res_diff, m4, 1e-9);
}

TEST(KernelGenerator, multi_result_kernel_reuse_kernel) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  matrix_cl<double> sum1_cl(3, 3);
  matrix_cl<double> diff1_cl(3, 3);
  matrix_cl<double> sum2_cl(3, 3);
  matrix_cl<double> diff2_cl(3, 3);

  auto sum = m1_cl + m2_cl;
  stan::math::results(sum1_cl, diff1_cl)
      = stan::math::expressions(sum, m1_cl - m2_cl);
  stan::math::results(sum2_cl, diff2_cl)
      = stan::math::expressions(sum, m1_cl - m2_cl);

  MatrixXd res_sum = stan::math::from_matrix_cl(sum2_cl);
  MatrixXd res_diff = stan::math::from_matrix_cl(diff2_cl);

  MatrixXd correct_sum = m1 + m2;
  MatrixXd correct_diff = m1 - m2;

  EXPECT_MATRIX_NEAR(res_sum, correct_sum, 1e-9);
  EXPECT_MATRIX_NEAR(res_diff, correct_diff, 1e-9);
}

#endif
