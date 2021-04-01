#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto unit_vector_constrain_functor
    = [](const auto& a) { return stan::math::unit_vector_constrain(a); };
auto unit_vector_constrain_functor2 = [](const auto& a) {
  using T_lp = stan::return_type_t<decltype(a)>;
  T_lp lp(4);
  if (!stan::is_constant<T_lp>::value) {
    stan::math::adjoint_of(lp) += 9;
  }
  return stan::math::unit_vector_constrain(a, lp);
};
auto unit_vector_constrain_functor3 = [](const auto& a) {
  using T_lp = stan::return_type_t<decltype(a)>;
  T_lp lp(4);
  stan::math::eval(stan::math::unit_vector_constrain(a, lp));
  return lp;
};

TEST(OpenCLUnitVectorConstrain, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;

  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor,
                                                a);
  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor2,
                                                a);
  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor3,
                                                a);
}

TEST(OpenCLUnitVectorConstrain, prim_rev_size_1) {
  int N = 1;

  Eigen::VectorXd a(N);
  a << 123;

  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor,
                                                a);
  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor2,
                                                a);
  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor3,
                                                a);
}

TEST(OpenCLUnitVectorConstrain, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);

  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor,
                                                a);
  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor2,
                                                a);
  stan::math::test::compare_cpu_opencl_prim_rev(unit_vector_constrain_functor3,
                                                a);
}

#endif
