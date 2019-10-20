#include <stan/math/rev/scal.hpp>
#include <stan/math/prim/scal.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

TEST(requires, var_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_t, var>::unary();
}
TEST(requires, var_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_var_t, var>::not_unary();
}
TEST(requires, var_all_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_var_t, var, var, var>::all();
}
TEST(requires, var_all_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_var_t, var, var>::all_not();
}
TEST(requires, var_any_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_var_t, var, var>::any();
}
TEST(requires, var_any_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_var_t, var, var>::any_not();
}

TEST(requires, var_or_fvar_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_or_fvar_t, var, var>::unary();
}
TEST(requires, var_or_fvar_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_var_or_fvar_t, var, var>::not_unary();
}
TEST(requires, var_or_fvar_all_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_var_or_fvar_t, var, var>::all();
}
TEST(requires, var_or_fvar_all_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_var_or_fvar_t, var,
                       var>::all_not();
}
TEST(requires, var_or_fvar_any_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_var_or_fvar_t, var, var>::any();
}
TEST(requires, var_or_fvar_any_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_var_or_fvar_t, var,
                       var>::any_not();
}

TEST(requires, stan_scalar_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_stan_scalar_t, var, var>::unary();
}
TEST(requires, stan_scalar_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_stan_scalar_t, var, var>::not_unary();
}
TEST(requires, stan_scalar_all_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_stan_scalar_t, var, var>::all();
}
TEST(requires, stan_scalar_all_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_stan_scalar_t, var,
                       var>::all_not();
}
TEST(requires, stan_scalar_any_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_stan_scalar_t, var, var>::any();
}
TEST(requires, stan_scalar_any_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_stan_scalar_t, var,
                       var>::any_not();
}
