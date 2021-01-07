#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto fabs_functor = [](const auto& a) { return stan::math::fabs(a); };

TEST(OpenCL_fabs, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -8, 2.7, 8, -2.6, -2, 1, 1.3, 3;
  stan::math::test::compare_cpu_opencl_prim_rev(fabs_functor, a);

  Eigen::VectorXd b_nan(3);
  b_nan << stan::math::NOT_A_NUMBER, stan::math::NOT_A_NUMBER,
      stan::math::NOT_A_NUMBER;
  stan::math::test::compare_cpu_opencl_prim_rev(fabs_functor, b_nan);
}

TEST(OpenCL_fabs, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(fabs_functor, a);
}

TEST(OpenCL_fabs, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(fabs_functor, a);
}

#endif
