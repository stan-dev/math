#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto plus_functor = [](const auto& a) { return stan::math::plus(a); };

TEST(OpenCLPlus, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(plus_functor, a);
}

#endif
