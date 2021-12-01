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

const char opencl_code_plus_minus[] = STRINGIFY(
    double pm_res1 = pm_in1 + pm_in2; double pm_res2 = pm_in1 - pm_in2;);

TEST(KernelGenerator, opencl_code_test) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::matrix_cl;
  VectorXd m1(6, 1);
  m1 << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;
  VectorXd m2(6, 1);
  m2 << 1.5, 1.6, 1.1, 1.2, 1.3, 1.4;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto code = stan::math::opencl_code<opencl_code_plus_minus>(
      std::make_tuple("pm_in1", "pm_in2"), m1_cl, m2_cl);

  matrix_cl<double> res1_cl;
  matrix_cl<double> res2_cl;

  stan::math::results(res1_cl, res2_cl) = stan::math::expressions(
      code.output<double>("pm_res1"), code.output<double>("pm_res2"));

  MatrixXd res1 = stan::math::from_matrix_cl(res1_cl);
  MatrixXd res2 = stan::math::from_matrix_cl(res2_cl);
  EXPECT_MATRIX_NEAR(m1 + m2, res1, 1e-9);
  EXPECT_MATRIX_NEAR(m1 - m2, res2, 1e-9);
}

const char opencl_code_constant[] = STRINGIFY(int res = 123450;);

TEST(KernelGenerator, opencl_code_no_inputs_test) {
  using Eigen::MatrixXi;
  using Eigen::VectorXi;
  using stan::math::matrix_cl;
  VectorXi m(6, 1);
  m << 1, 2, 3, 4, 5, 6;
  VectorXi correct = m.array() + 123450;

  matrix_cl<int> m_cl(m);

  auto code = stan::math::opencl_code<opencl_code_constant>(std::make_tuple());

  matrix_cl<int> res_cl = m_cl + code.output<int>("res");

  MatrixXi res = stan::math::from_matrix_cl(res_cl);
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

#endif
