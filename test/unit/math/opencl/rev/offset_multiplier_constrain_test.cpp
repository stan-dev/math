#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto offset_multiplier_constrain_functor
    = [](const auto& a, const auto& b, const auto& c) {
        return stan::math::offset_multiplier_constrain(a, b, c);
      };
auto offset_multiplier_constrain_functor2
    = [](const auto& a, const auto& b, const auto& c) {
        using T_lp = stan::return_type_t<decltype(a), decltype(b), decltype(c)>;
        T_lp lp(4);
        if (!stan::is_constant<T_lp>::value) {
          stan::math::adjoint_of(lp) += 9;
        }
        return stan::math::offset_multiplier_constrain(a, b, c, lp);
      };
auto offset_multiplier_constrain_functor3
    = [](const auto& a, const auto& b, const auto& c) {
        using T_lp = stan::return_type_t<decltype(a), decltype(b), decltype(c)>;
        T_lp lp(4);
        stan::math::eval(stan::math::offset_multiplier_constrain(a, b, c, lp));
        return lp;
      };

TEST(OpenCLOffsetMultiplierConstrain, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::VectorXd b(8);
  b << 0.5, 1, 1.5, 2.2, 3, 4, 3.4, 0.8;
  Eigen::VectorXd c(8);
  c << 3, 43.4, 2.2, 0.8, 0.5, 1, 1.5, 3;
  double b_scal = 0.9;
  double c_scal = 12.4;

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b, c);

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b_scal, c);

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b, c_scal);
}

TEST(OpenCLOffsetMultiplierConstrain, prim_rev_size_0) {
  int N = 0;

  Eigen::VectorXd a(N);
  Eigen::VectorXd b(N);
  Eigen::VectorXd c(N);
  double b_scal = 0.9;
  double c_scal = 12.4;

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b, c);

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b_scal, c);

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b, c_scal);
}

TEST(OpenCLOffsetMultiplierConstrain, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N).array().abs();
  Eigen::VectorXd c = Eigen::VectorXd::Random(N).array().abs();
  double b_scal = 0.9;
  double c_scal = 12.4;

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b, c);

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b_scal, c);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b_scal, c);

  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor, a, b, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor2, a, b, c_scal);
  stan::math::test::compare_cpu_opencl_prim_rev(
      offset_multiplier_constrain_functor3, a, b, c_scal);
}

#endif
