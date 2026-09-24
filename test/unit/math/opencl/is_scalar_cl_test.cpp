#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST(type_trait, is_scalar_cl) {
  using stan::is_scalar_cl;
  using stan::math::matrix_cl;
  using stan::math::var;
  using stan::math::opencl::ScalarCl;
  EXPECT_TRUE((is_scalar_cl<ScalarCl<double>>::value));
  EXPECT_TRUE((is_scalar_cl<ScalarCl<var>>::value));
  EXPECT_TRUE((is_scalar_cl<const ScalarCl<double>&>::value));
  EXPECT_TRUE((is_scalar_cl<ScalarCl<var>&&>::value));
  EXPECT_FALSE((is_scalar_cl<double>::value));
  EXPECT_FALSE((is_scalar_cl<var>::value));
  EXPECT_FALSE((is_scalar_cl<matrix_cl<double>>::value));
}

TEST(type_trait, scalar_cl_value_and_scalar_type) {
  using stan::scalar_type_t;
  using stan::value_type_t;
  using stan::math::var;
  using stan::math::opencl::ScalarCl;
  EXPECT_TRUE((std::is_same<scalar_type_t<ScalarCl<double>>, double>::value));
  EXPECT_TRUE((std::is_same<value_type_t<ScalarCl<double>>, double>::value));
}

TEST(type_trait, prim_and_rev_scalar_cl) {
  using stan::is_prim_scalar_cl;
  using stan::is_rev_scalar_cl;
  using stan::math::var;
  using stan::math::opencl::ScalarCl;
  EXPECT_TRUE((is_prim_scalar_cl<ScalarCl<double>>::value));
  EXPECT_TRUE((is_prim_scalar_cl<const ScalarCl<double>&>::value));
  EXPECT_FALSE((is_prim_scalar_cl<ScalarCl<var>>::value));
  EXPECT_FALSE((is_prim_scalar_cl<double>::value));
  EXPECT_TRUE((is_rev_scalar_cl<ScalarCl<var>>::value));
  EXPECT_TRUE((is_rev_scalar_cl<ScalarCl<var>&>::value));
  EXPECT_FALSE((is_rev_scalar_cl<ScalarCl<double>>::value));
  EXPECT_FALSE((is_rev_scalar_cl<var>::value));
}

TEST(type_trait, scalar_cl_is_not_a_host_scalar) {
  using stan::is_stan_scalar;
  using stan::is_var;
  using stan::math::var;
  using stan::math::opencl::ScalarCl;
  EXPECT_FALSE((is_stan_scalar<ScalarCl<double>>::value));
  EXPECT_FALSE((std::is_arithmetic<ScalarCl<double>>::value));
}

TEST(type_trait, scalar_cl_is_autodiff) {
  using stan::is_autodiff_v;
  using stan::math::var;
  using stan::math::opencl::ScalarCl;
  EXPECT_FALSE((is_autodiff_v<ScalarCl<double>>));
}
#endif
