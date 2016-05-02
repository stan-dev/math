#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_EQ_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_EQ_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <stan/math/fwd/mat.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>

static inline void expect_val_deriv_eq(double exp_var, double test_var) {
  expect_val_eq(exp_var, test_var);
}

template <typename T>
static inline void expect_val_deriv_eq(T exp_var, T test_var) {
  expect_val_deriv_eq(exp_var.val(), test_var.val());
  expect_val_deriv_eq(exp_var.d_, test_var.d_);
}

#endif
