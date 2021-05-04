#include <stan/math/rev/meta/is_var.hpp>  // just this is bad

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
#include <string>

TEST(requires_rev_scal, var_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_t, var>::unary();
}
TEST(requires_rev_scal, var_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_var_t, var>::not_unary();
}
TEST(requires_rev_scal, var_all_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_var_t, var, var, var>::all();
}
TEST(requires_rev_scal, var_all_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_var_t, var, var>::all_not();
}
TEST(requires_rev_scal, var_any_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_var_t, var, var>::any();
}
TEST(requires_rev_scal, var_any_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_var_t, var, var>::any_not();
}

TEST(requires_rev_scal, autodiff_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_autodiff_t, var, var>::unary();
}
TEST(requires_rev_scal, autodiff_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_autodiff_t, var, var>::not_unary();
}
TEST(requires_rev_scal, autodiff_all_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_autodiff_t, var, var>::all();
}
TEST(requires_rev_scal, autodiff_all_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_autodiff_t, var, var>::all_not();
}
TEST(requires_rev_scal, autodiff_any_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_autodiff_t, var, var>::any();
}
TEST(requires_rev_scal, autodiff_any_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_autodiff_t, var, var>::any_not();
}

TEST(requires_rev_scal, stan_scalar_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_stan_scalar_t, var, var>::unary();
}
TEST(requires_rev_scal, stan_scalar_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_stan_scalar_t, var, var>::not_unary();
}
TEST(requires_rev_scal, stan_scalar_all_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_stan_scalar_t, var, var>::all();
}
TEST(requires_rev_scal, stan_scalar_all_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_stan_scalar_t, var,
                       var>::all_not();
}
TEST(requires_rev_scal, stan_scalar_any_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_stan_scalar_t, var, var>::any();
}
TEST(requires_rev_scal, stan_scalar_any_not_test) {
  using stan::math::var;
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_stan_scalar_t, var,
                       var>::any_not();
}

template <typename... Types>
using require_test_var_return_t
    = stan::require_return_type_t<stan::is_var, Types...>;

template <typename... Types>
using require_test_arithmetic_return_t
    = stan::require_return_type_t<std::is_arithmetic, Types...>;

template <typename... Types>
using require_test_not_var_return_t
    = stan::require_not_return_type_t<stan::is_var, Types...>;

template <typename... Types>
using require_test_not_arithmetic_return_t
    = stan::require_not_return_type_t<std::is_arithmetic, Types...>;

TEST(requires_rev_scal, return_type_t_test) {
  using stan::require_return_type_t;
  using stan::math::var;
  using stan::math::var_value;
  using stan::test::require_variadic_checker;

  EXPECT_TRUE((
      require_variadic_checker<require_test_var_return_t, var, double>::value));
  EXPECT_TRUE((require_variadic_checker<require_test_var_return_t, var, double,
                                        Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((require_variadic_checker<require_test_arithmetic_return_t, var,
                                         double>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_test_arithmetic_return_t, var, double,
                                Eigen::Matrix<var, -1, -1>>::value));

  EXPECT_FALSE((require_variadic_checker<require_test_not_var_return_t, var,
                                         double>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_test_not_var_return_t, var, double,
                                Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_TRUE((require_variadic_checker<require_test_not_arithmetic_return_t,
                                        var, double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_test_not_arithmetic_return_t, var,
                                double, Eigen::Matrix<var, -1, -1>>::value));
}
