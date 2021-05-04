#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto minus_functor = [](const auto& a) { return stan::math::minus(a); };

TEST(OpenCLMinus, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(minus_functor, a);
}

#endif
