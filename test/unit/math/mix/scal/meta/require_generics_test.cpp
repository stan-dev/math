#include <stan/math/fwd/scal.hpp>
#include <stan/math/rev/scal.hpp>
#include <stan/math/prim/scal.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>


TEST(requires, var_or_fvar_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_var_or_fvar, var, fvar<double>>::unary();
}
TEST(requires, var_or_fvar_not_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_not_var_or_fvar, var, fvar<double>>::not_unary();
}
TEST(requires, var_or_fvar_all_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_all_var_or_fvar, var, fvar<double>>::all();
}
TEST(requires, var_or_fvar_all_not_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_all_not_var_or_fvar, var, fvar<double>>::all_not();
}
TEST(requires, var_or_fvar_any_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_any_var_or_fvar, var, fvar<double>>::any();
}
TEST(requires, var_or_fvar_any_not_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_any_not_var_or_fvar, var, fvar<double>>::any_not();
}

TEST(requires, stan_scalar_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_stan_scalar, var, fvar<double>>::unary();
}
TEST(requires, stan_scalar_not_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_not_stan_scalar, var, fvar<double>>::not_unary();
}
TEST(requires, stan_scalar_all_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_all_stan_scalar, var, fvar<double>>::all();
}
TEST(requires, stan_scalar_all_not_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_all_not_stan_scalar, var, fvar<double>>::all_not();
}
TEST(requires, stan_scalar_any_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_any_stan_scalar, var, fvar<double>>::any();
}
TEST(requires, stan_scalar_any_not_mix_test) {
  using stan::test::require_autodiff_checker;
  using stan::math::fvar; using stan::math::var;
  require_autodiff_checker<stan::require_any_not_stan_scalar, var, fvar<double>>::any_not();
}
