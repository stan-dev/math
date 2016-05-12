#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_VAL_DERIV_EQ_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_VAL_DERIV_EQ_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_val_deriv_eq.hpp>
#include <gtest/gtest.h>

template <typename V>
static inline void expect_val_deriv_eq(V exp_var, V base_exp_var, 
                                       V test_var, V base_test_var) {
  expect_val_deriv_eq(exp_var.val(), base_exp_var.val(),
                      test_var.val(), base_test_var.val());
  expect_val_deriv_eq(exp_var.d_, base_exp_var.d_, 
                      test_var.d_, base_test_var.d_);
}

#endif
