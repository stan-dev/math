#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto log_softmax_functor
    = [](const auto& a) { return stan::math::log_softmax(a); };

TEST(OpenCLLogSoftmax, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1 + std::numeric_limits<double>::epsilon(), 1.5, 3, 3.4,
      4;
  stan::math::test::compare_cpu_opencl_prim_rev(log_softmax_functor, a);
}

TEST(OpenCLLogSoftmax, prim_rev_size_1) {
  int N = 1;

  Eigen::VectorXd a(N);
  a << 23;
  stan::math::test::compare_cpu_opencl_prim_rev(log_softmax_functor, a);
}

TEST(OpenCLLogSoftmax, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(log_softmax_functor, a);
}

#endif
