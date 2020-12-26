#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto mean_functor = [](const auto& a) { return stan::math::mean(a); };

TEST(OpenCLMean, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(mean_functor, a);

  Eigen::MatrixXd b(1, 1);
  b << 2.5;
  stan::math::test::compare_cpu_opencl_prim_rev(mean_functor, b);
}

TEST(OpenCLMean, exceptions) {
  using stan::math::matrix_cl;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, -1> m0(0, 0);
  matrix_cl<double> m0_cl = stan::math::to_matrix_cl(m0);
  var_value<matrix_cl<double>> m0_var_cl = m0_cl;

  EXPECT_THROW(stan::math::mean(m0_cl), std::invalid_argument);
  EXPECT_THROW(stan::math::mean(m0_var_cl), std::invalid_argument);
}

#endif
