#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/err.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(KernelGenerator, check_cl_error) {
  Eigen::MatrixXd m1(1, 3);
  m1 << 1, 2, 0.1;
  stan::math::matrix_cl<double> m1_cl(m1);
  stan::math::matrix_cl<double> m2_cl(3, 2);
  EXPECT_THROW(
      stan::math::check_cl("test", "m1", m1_cl, "positive") = m2_cl > 0,
      std::invalid_argument);
}

TEST(KernelGenerator, check_cl_scalar_check_only_pass) {
  EXPECT_NO_THROW(stan::math::check_cl("test", "m1", 3, "positive") = 3 > 0);
}

TEST(KernelGenerator, check_cl_scalar_check_only_throw) {
  EXPECT_THROW(stan::math::check_cl("test", "m1", 0, "positive") = 0 > 0,
               std::domain_error);
}

TEST(KernelGenerator, check_cl_check_only_pass) {
  std::string kernel_filename = "check_cl_positive.cl";
  Eigen::MatrixXd m1(1, 3);
  m1 << 1, 2, 0.1;
  stan::math::matrix_cl<double> m1_cl(m1);
  auto res = stan::math::check_cl("test", "m1", m1_cl, "positive");
  auto expr = m1_cl > 0;
  std::string kernel_src = expr.get_kernel_source_for_evaluating_into(res);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  EXPECT_NO_THROW(res = expr);
}

TEST(KernelGenerator, check_cl_check_only_throw) {
  Eigen::MatrixXd m1(1, 3);
  m1 << 1, -1, 0;
  stan::math::matrix_cl<double> m1_cl(m1);
  EXPECT_THROW(
      stan::math::check_cl("test", "m1", m1_cl, "positive") = m1_cl > 0,
      std::domain_error);
}

TEST(KernelGenerator, check_cl_in_kernel_pass) {
  Eigen::MatrixXd m1(1, 3);
  m1 << 1, 2, 0.1;
  stan::math::matrix_cl<double> m1_cl(m1);
  stan::math::matrix_cl<double> m2_cl;
  EXPECT_NO_THROW(
      stan::math::results(m2_cl,
                          stan::math::check_cl("test", "m1", m1_cl, "positive"),
                          stan::math::check_cl("test", "scal", 3, "positive"))
      = stan::math::expressions(m1_cl * 2 + 1, m1_cl > 0, 3 > 0));
}

TEST(KernelGenerator, check_cl_in_kernel_throw) {
  Eigen::MatrixXd m1(1, 3);
  m1 << 1, -1, 0;
  stan::math::matrix_cl<double> m1_cl(m1);
  stan::math::matrix_cl<double> m2_cl;
  EXPECT_THROW(stan::math::results(
                   m2_cl, stan::math::check_cl("test", "m1", m1_cl, "positive"))
               = stan::math::expressions(m1_cl * 2 + 1, m1_cl > 0),
               std::domain_error);
}
#endif
