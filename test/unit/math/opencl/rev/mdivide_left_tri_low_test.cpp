#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto mdivide_left_tri_low_functor1 = [](const auto& a, const auto& b) {
  return stan::math::mdivide_left_tri_low(a, b);
};

auto mdivide_left_tri_low_functor2
    = [](const auto& a) { return stan::math::mdivide_left_tri_low(a); };

TEST(OpenCL_mdivide_left_tri_low, mdivide_left_tri_low_small) {
  Eigen::MatrixXd in1(3, 3);
  in1 << 0.5, 3.4, 5.2, 7.5, 1, 2, 3, 4, 5;
  Eigen::VectorXd in2(3);
  in2 << 3.3, 0.9, 6.7;
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_left_tri_low_functor1,
                                                in1, in2);
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_left_tri_low_functor2,
                                                in1);
}

TEST(OpenCL_mdivide_left_tri_low, zero) {
  Eigen::MatrixXd in1;
  Eigen::VectorXd in2;
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_left_tri_low_functor1,
                                                in1, in2);
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_left_tri_low_functor2,
                                                in1);
}

TEST(OpenCL_mdivide_left_tri_low, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_left_tri_low_functor1,
                                                a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_left_tri_low_functor2,
                                                a);
}

#endif
