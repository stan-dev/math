#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_FVAR_VAR_EQ_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_FVAR_VAR_EQ_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>

static inline void expect_fvar_var_eq(stan::math::var exp_var,
                               stan::math::var base_exp_var,
                                 stan::math::var test_var,
                                   stan::math::var base_test_var) {
  EXPECT_FLOAT_EQ(exp_var.val(), test_var.val());
  AVEC exp_y = createAVEC(base_exp_var);
  VEC exp_g;
  exp_var.grad(exp_y, exp_g);
  AVEC test_y = createAVEC(base_test_var);
  VEC test_g;
  test_var.grad(test_y, test_g);
  EXPECT_FLOAT_EQ(exp_g[0], test_g[0]);
}

template <typename V>
static inline void expect_fvar_var_eq(V exp_var, V base_exp_var, 
                               V test_var, V base_test_var) {
  expect_fvar_var_eq(exp_var.val(), base_exp_var.val(),
              test_var.val(), base_test_var.val());
  expect_fvar_var_eq(exp_var.d_, base_exp_var.d_, 
              test_var.d_, base_test_var.d_);
}

#endif
