#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto acoshh_functor = [](const auto& a) { return stan::math::acosh(a); };

TEST(OpenCLacoshh, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << 2.2, 1.8, 1, 1 + std::numeric_limits<double>::epsilon(), 1.2, 1.3, 1.34,
      1.4;
  stan::math::test::compare_cpu_opencl_prim_rev(acoshh_functor, a);
}

TEST(OpenCLacoshh, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(acoshh_functor, a);
}

TEST(OpenCLacoshh, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a
      = Eigen::MatrixXd::Constant(N, N, 2.0) + Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(acoshh_functor, a);
}

#endif
