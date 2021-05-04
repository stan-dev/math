#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto atanh_functor = [](const auto& a) { return stan::math::atanh(a); };

TEST(OpenCLAtanh, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -0.98, -0.8, 0.5, 0.6 + std::numeric_limits<double>::epsilon(), 0.5,
      -0.03, -0.34, 0.44;
  stan::math::test::compare_cpu_opencl_prim_rev(atanh_functor, a);
}

TEST(OpenCLAtanh, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(atanh_functor, a);
}

TEST(OpenCLAtanh, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(atanh_functor, a);
}

#endif
