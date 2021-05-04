#include <stan/math/fwd.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>

TEST(requires_fwd_scal, fvar_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_fvar_t, fvar<double>>::unary();
}
TEST(requires_fwd_scal, fvar_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_fvar_t, fvar<double>>::not_unary();
}
TEST(requires_fwd_scal, fvar_all_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_fvar_t, fvar<double>, fvar<double>,
                       fvar<double>>::all();
}
TEST(requires_fwd_scal, fvar_all_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_fvar_t, fvar<double>,
                       fvar<double>>::all_not();
}
TEST(requires_fwd_scal, fvar_any_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_fvar_t, fvar<double>,
                       fvar<double>>::any();
}
TEST(requires_fwd_scal, fvar_any_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_fvar_t, fvar<double>,
                       fvar<double>>::any_not();
}

TEST(requires_fwd_scal, autodiff_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_autodiff_t, fvar<double>,
                       fvar<double>>::unary();
}
TEST(requires_fwd_scal, autodiff_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_autodiff_t, fvar<double>,
                       fvar<double>>::not_unary();
}
TEST(requires_fwd_scal, autodiff_all_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_autodiff_t, fvar<double>,
                       fvar<double>>::all();
}
TEST(requires_fwd_scal, autodiff_all_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_autodiff_t, fvar<double>,
                       fvar<double>>::all_not();
}
TEST(requires_fwd_scal, autodiff_any_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_autodiff_t, fvar<double>,
                       fvar<double>>::any();
}
TEST(requires_fwd_scal, autodiff_any_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_autodiff_t, fvar<double>,
                       fvar<double>>::any_not();
}

TEST(requires_fwd_scal, stan_scalar_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_stan_scalar_t, fvar<double>,
                       fvar<double>>::unary();
}
TEST(requires_fwd_scal, stan_scalar_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_stan_scalar_t, fvar<double>,
                       fvar<double>>::not_unary();
}
TEST(requires_fwd_scal, stan_scalar_all_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_stan_scalar_t, fvar<double>,
                       fvar<double>>::all();
}
TEST(requires_fwd_scal, stan_scalar_all_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_stan_scalar_t, fvar<double>,
                       fvar<double>>::all_not();
}
TEST(requires_fwd_scal, stan_scalar_any_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_stan_scalar_t, fvar<double>,
                       fvar<double>>::any();
}
TEST(requires_fwd_scal, stan_scalar_any_not_test) {
  using stan::math::fvar;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_stan_scalar_t, fvar<double>,
                       fvar<double>>::any_not();
}
