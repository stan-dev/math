#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto rows_dot_self_functor
    = [](const auto& a) { return stan::math::rows_dot_self(a); };

TEST(OpenCLRowsDotSelf, prim_rev_small_vector) {
  Eigen::MatrixXd a(2, 4);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(rows_dot_self_functor, a);
}

TEST(OpenCLRowsDotSelf, prim_rev_size_0) {
  Eigen::MatrixXd a(0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(rows_dot_self_functor, a);
}

TEST(OpenCLRowsDotSelf, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(rows_dot_self_functor, a);
}

#endif
