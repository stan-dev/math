#include <stan/math/fwd.hpp>
#include <stan/math/rev/scal.hpp>
#include <stan/math/prim/scal.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

TEST(requires_mix_scal, var_or_fvar_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_or_fvar_t, var, fvar<double>>::unary();
}
TEST(requires_mix_scal, var_or_fvar_not_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_var_or_fvar_t, var,
                       fvar<double>>::not_unary();
}
TEST(requires_mix_scal, var_or_fvar_all_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_var_or_fvar_t, var,
                       fvar<double>>::all();
}
TEST(requires_mix_scal, var_or_fvar_all_not_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_var_or_fvar_t, var,
                       fvar<double>>::all_not();
}
TEST(requires_mix_scal, var_or_fvar_any_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_var_or_fvar_t, var,
                       fvar<double>>::any();
}
TEST(requires_mix_scal, var_or_fvar_any_not_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_var_or_fvar_t, var,
                       fvar<double>>::any_not();
}

TEST(requires_mix_scal, stan_scalar_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_stan_scalar_t, var, fvar<double>>::unary();
}
TEST(requires_mix_scal, stan_scalar_not_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_stan_scalar_t, var,
                       fvar<double>>::not_unary();
}
TEST(requires_mix_scal, stan_scalar_all_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_stan_scalar_t, var,
                       fvar<double>>::all();
}
TEST(requires_mix_scal, stan_scalar_all_not_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_stan_scalar_t, var,
                       fvar<double>>::all_not();
}
TEST(requires_mix_scal, stan_scalar_any_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_stan_scalar_t, var,
                       fvar<double>>::any();
}
TEST(requires_mix_scal, stan_scalar_any_not_mix_test) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_stan_scalar_t, var,
                       fvar<double>>::any_not();
}
