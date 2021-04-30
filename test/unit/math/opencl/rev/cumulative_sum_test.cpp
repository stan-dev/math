#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto cumulative_sum_functor
    = [](const auto& a) { return stan::math::cumulative_sum(a); };

TEST(OpenCLPrim, cumulative_sum_small) {
  Eigen::VectorXd d(11);
  d << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;
  stan::math::test::compare_cpu_opencl_prim(cumulative_sum_functor, d);
}

TEST(OpenCLPrim, cumulative_sum_sizee0) {
  Eigen::VectorXd d(0);
  stan::math::test::compare_cpu_opencl_prim(cumulative_sum_functor, d);
}

TEST(OpenCLPrim, cumulative_sum_sizee1) {
  Eigen::VectorXd d(1);
  d << 12;
  stan::math::test::compare_cpu_opencl_prim(cumulative_sum_functor, d);
}

TEST(OpenCLPrim, cumulative_sum_large) {
  Eigen::VectorXd d = Eigen::VectorXd::Random(500000);
  stan::math::test::compare_cpu_opencl_prim(cumulative_sum_functor, d);
}

#endif
