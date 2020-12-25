#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto lbeta_functor
    = [](const auto& a, const auto& b) { return stan::math::lbeta(a, b); };

TEST(OpenCL_lbeta, value_small) {
  Eigen::VectorXd in1(4);
  in1 << 0.5, 3.4, 5.2, 7.5;
  Eigen::VectorXd in2(4);
  in2 << 3.3, 0.9, 6.7, 1.8;
  stan::math::test::compare_cpu_opencl_prim_rev(lbeta_functor, in1, in2);
}

#endif
