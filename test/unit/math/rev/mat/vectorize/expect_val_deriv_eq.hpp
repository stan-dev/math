#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VAL_DERIV_EQ_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VAL_DERIV_EQ_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>

static inline void expect_val_deriv_eq(stan::math::var exp_var,
                                       stan::math::var base_exp_var,
                                       stan::math::var test_var,
                                       stan::math::var base_test_var) {
  using stan::math::set_zero_all_adjoints;

  expect_val_eq(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  expect_val_eq(exp_g[0], test_g[0]);
}

#endif
