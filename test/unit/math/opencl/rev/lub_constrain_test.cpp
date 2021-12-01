#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto lub_constrain_functor = [](const auto& a, const auto& b, const auto& c) {
  return stan::math::lub_constrain(a, b, c);
};
auto lub_constrain_functor2 = [](const auto& a, const auto& b, const auto& c) {
  using T_lp = stan::return_type_t<decltype(a), decltype(b), decltype(c)>;
  T_lp lp(4);
  if (!stan::is_constant<T_lp>::value) {
    stan::math::adjoint_of(lp) += 9;
  }
  return stan::math::lub_constrain(a, b, c, lp);
};
auto lub_constrain_functor3 = [](const auto& a, const auto& b, const auto& c) {
  using T_lp = stan::return_type_t<decltype(a), decltype(b), decltype(c)>;
  T_lp lp(4);
  stan::math::eval(stan::math::lub_constrain(a, b, c, lp));
  return lp;
};

TEST(OpenCLLubConstrain, prim_rev_values_small) {
  Eigen::VectorXd a(7);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4;
  Eigen::VectorXd b(7);
  b << -2.2, -0.8, 1, 4, -INFINITY, -INFINITY, 3;
  Eigen::VectorXd c(7);
  c << -2.1, 0.8, 4, 7, INFINITY, 2, INFINITY;
  double b_scal = -8;
  double c_scal = 8;

  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c_scal);
}

TEST(OpenCLLubConstrain, prim_rev_size_0) {
  int N = 0;

  Eigen::RowVectorXd a(N);
  Eigen::RowVectorXd b(N);
  Eigen::RowVectorXd c(N);
  double b_scal = -8;
  double c_scal = 8;

  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c_scal);
}

TEST(OpenCLLubConstrain, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, N);
  Eigen::MatrixXd c = b.array() + Eigen::ArrayXXd::Random(N, N) + 1.0;
  double b_scal = -1.5;
  double c_scal = 1.5;

  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor2, a,
                                                b_scal, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a, b,
                                                c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(lub_constrain_functor3, a,
                                                b_scal, c_scal);
}

#endif
