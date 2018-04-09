#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_BINARY_VAL_DERIV_EQ_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_BINARY_VAL_DERIV_EQ_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <gtest/gtest.h>

static inline void expect_binary_val_deriv_eq(stan::math::var exp_var,
                                              stan::math::var base_exp_var1,
                                              stan::math::var base_exp_var2,
                                              stan::math::var test_var,
                                              stan::math::var base_test_var1,
                                              stan::math::var base_test_var2) {
  using stan::math::set_zero_all_adjoints;
  expect_val_eq(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var1, base_exp_var2);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var1, base_test_var2);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  for (size_t i = 0; i < exp_g.size(); ++i)
    expect_val_eq(exp_g[i], test_g[i]);
}

static inline void expect_binary_val_deriv_eq(stan::math::var exp_var,
                                              stan::math::var base_exp_var1,
                                              double base_exp_var2,
                                              stan::math::var test_var,
                                              stan::math::var base_test_var1,
                                              double base_test_var2) {
  using stan::math::set_zero_all_adjoints;
  expect_val_eq(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var1);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var1);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  expect_val_eq(exp_g[0], test_g[0]);
}

static inline void expect_binary_val_deriv_eq(stan::math::var exp_var,
                                              double base_exp_var1,
                                              stan::math::var base_exp_var2,
                                              stan::math::var test_var,
                                              double base_test_var1,
                                              stan::math::var base_test_var2) {
  using stan::math::set_zero_all_adjoints;
  expect_val_eq(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var2);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var2);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  expect_val_eq(exp_g[0], test_g[0]);
}

static inline void expect_binary_val_deriv_eq(stan::math::var exp_var,
                                              stan::math::var base_exp_var1,
                                              int base_exp_var2,
                                              stan::math::var test_var,
                                              stan::math::var base_test_var1,
                                              int base_test_var2) {
  using stan::math::set_zero_all_adjoints;
  expect_val_eq(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var1);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var1);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  expect_val_eq(exp_g[0], test_g[0]);
}

static inline void expect_binary_val_deriv_eq(stan::math::var exp_var,
                                              int base_exp_var1,
                                              stan::math::var base_exp_var2,
                                              stan::math::var test_var,
                                              int base_test_var1,
                                              stan::math::var base_test_var2) {
  using stan::math::set_zero_all_adjoints;
  expect_val_eq(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var2);
  VEC exp_g;
  set_zero_all_adjoints();
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var2);
  VEC test_g;
  set_zero_all_adjoints();
  test_var.grad(test_y, test_g);
  expect_val_eq(exp_g[0], test_g[0]);
}
#endif
