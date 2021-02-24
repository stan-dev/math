#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto distance_functor
    = [](const auto& a, const auto& b) { return stan::math::distance(a, b); };

TEST(OpenCLMatrix_distance, prim_rev_values_small) {
  int N = 6;

  Eigen::VectorXd a(N);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::VectorXd b(N);
  b << 12, 11, 10, 9, 8, 7;
  stan::math::test::compare_cpu_opencl_prim_rev(distance_functor, a, b);
}

TEST(OpenCLMatrix_distance, prim_rev_values_M_0) {
  int N = 0;

  Eigen::VectorXd a(N);
  Eigen::VectorXd b(N);
  stan::math::test::compare_cpu_opencl_prim_rev(distance_functor, a, b);
}

TEST(OpenCLMatrix_distance, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(distance_functor, a, b);
}

#endif
