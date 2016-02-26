#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_EQ_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_EQ_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <gtest/gtest.h>

static inline void expect_val_deriv_eq(double exp_var, double test_var) {
  EXPECT_FLOAT_EQ(exp_var, test_var);
}

template <typename T>
static inline void expect_val_deriv_eq(T exp_var, T test_var) {
  expect_val_deriv_eq(exp_var.val(), test_var.val());
  expect_val_deriv_eq(exp_var.d_, test_var.d_);
}

#endif
