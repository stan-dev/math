#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto reverse_functor = [](const auto& a) { return stan::math::reverse(a); };

TEST(OpenCLReverse, error) {
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 76, 8, 9;
  stan::math::matrix_cl<double> a_cl(a);
  EXPECT_THROW(stan::math::reverse(a_cl), std::invalid_argument);
}

TEST(OpenCLReverse, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1 + std::numeric_limits<double>::epsilon(), 1.5, 3, 3.4,
      4;
  stan::math::test::compare_cpu_opencl_prim_rev(reverse_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(
      reverse_functor, stan::math::eval(stan::math::transpose(a)));
}

TEST(OpenCLReverse, prim_rev_size_0) {
  int N = 0;

  Eigen::VectorXd a(N);
  stan::math::test::compare_cpu_opencl_prim_rev(reverse_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(
      reverse_functor, stan::math::eval(stan::math::transpose(a)));
}

TEST(OpenCLReverse, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(reverse_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(
      reverse_functor, stan::math::eval(stan::math::transpose(a)));
}

#endif
