#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto lb_constrain_functor = [](const auto& a, const auto& b) {
  return stan::math::lb_constrain(a, b);
};
auto lb_constrain_functor2 = [](const auto& a, const auto& b) {
  using T_lp = stan::return_type_t<decltype(a), decltype(b)>;
  T_lp lp(4);
  if (!stan::is_constant<T_lp>::value) {
    stan::math::adjoint_of(lp) += 9;
  }
  return stan::math::lb_constrain(a, b, lp);
};
auto lb_constrain_functor3 = [](const auto& a, const auto& b) {
  using T_lp = stan::return_type_t<decltype(a), decltype(b)>;
  T_lp lp(4);
  stan::math::eval(stan::math::lb_constrain(a, b, lp));
  return lp;
};

TEST(OpenCLLbConstrain, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::VectorXd b(8);
  b << INFINITY, -INFINITY, 1, 2, 3, 4, INFINITY, -INFINITY;
  double c = 8;

  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor, a, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor2, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor2, a, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor3, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor3, a, c);
}

TEST(OpenCLLbConstrain, prim_rev_size_0) {
  int N = 0;

  Eigen::RowVectorXd a(N);
  Eigen::RowVectorXd b(N);
  double c = 5;

  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor, a, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor2, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor2, a, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor3, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor3, a, c);
}

TEST(OpenCLLbConstrain, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, N);
  double c = 5;

  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor, a, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor2, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor2, a, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor3, a, b);
  stan::math::test::compare_cpu_opencl_prim_rev(lb_constrain_functor3, a, c);
}

#endif
