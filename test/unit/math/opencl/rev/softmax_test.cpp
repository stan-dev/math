#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto softmax_functor = [](const auto& a) { return stan::math::softmax(a); };

TEST(OpenCLSoftmax, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(softmax_functor, a);
}

TEST(OpenCLSoftmax, prim_rev_size_0) {
  int N = 0;

  Eigen::VectorXd a(N);
  stan::math::test::compare_cpu_opencl_prim_rev(softmax_functor, a);
}

TEST(OpenCLSoftmax, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(softmax_functor, a);
}

#endif
