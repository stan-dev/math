#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto mdivide_right_tri_low_functor = [](const auto& a, const auto& b) {
  return stan::math::mdivide_right_tri_low(a, b);
};

TEST(OpenCL_mdivide_right_tri_low, mdivide_right_tri_low_small) {
  Eigen::MatrixXd in1(3, 3);
  in1 << 0.5, 3.4, 5.2, 7.5, 1, 2, 3, 4, 5;
  Eigen::RowVectorXd in2(3);
  in2 << 3.3, 0.9, 6.7;
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_right_tri_low_functor,
                                                in2, in1);
}

TEST(OpenCL_mdivide_right_tri_low, zero) {
  Eigen::MatrixXd in1;
  Eigen::RowVectorXd in2;
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_right_tri_low_functor,
                                                in2, in1);
}

TEST(OpenCL_mdivide_right_tri_low, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  Eigen::RowVectorXd b = Eigen::RowVectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(mdivide_right_tri_low_functor,
                                                b, a);
}

#endif
