#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_EQ_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_EQ_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>

static inline void expect_eq(stan::math::var exp_var,
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

#endif
