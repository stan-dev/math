#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VAL_DERIV_EQ_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VAL_DERIV_EQ_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/is_nan.hpp>

static inline void expect_val_deriv_eq(stan::math::var exp_var,
                                       stan::math::var base_exp_var,
                                       stan::math::var test_var,
                                       stan::math::var base_test_var) {
  using stan::math::set_zero_all_adjoints;
  using stan::math::is_nan;

  if (!(is_nan(exp_var.val()) && is_nan(test_var.val())))
    EXPECT_FLOAT_EQ(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  if (is_nan(exp_g[0]) && is_nan(test_g[0])) return;
  EXPECT_FLOAT_EQ(exp_g[0], test_g[0]);
}

#endif
