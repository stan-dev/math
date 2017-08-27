#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_BINARY_VAL_DERIV_EQ_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_BINARY_VAL_DERIV_EQ_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>

template <typename FV>
static inline void expect_binary_val_deriv_eq(
    stan::math::fvar<FV> exp_var, stan::math::fvar<FV> base_exp_var1,
    stan::math::fvar<FV> base_exp_var2, stan::math::fvar<FV> test_var,
    stan::math::fvar<FV> base_test_var1,
    stan::math::fvar<FV> base_test_var2) {
  expect_binary_val_deriv_eq(exp_var.val(), base_exp_var1.val(),
                             base_exp_var2.val(), test_var.val(),
                             base_test_var1.val(), base_test_var2.val());
  expect_binary_val_deriv_eq(exp_var.d_, base_exp_var1.d_,
                             base_exp_var2.d_, test_var.d_,
                             base_test_var1.d_, base_test_var2.d_);
}

template <typename FV>
static inline void expect_binary_val_deriv_eq(
    stan::math::fvar<FV> exp_var, stan::math::fvar<FV> base_exp_var1,
    int base_exp_var2, stan::math::fvar<FV> test_var,
    stan::math::fvar<FV> base_test_var1, int base_test_var2) {
  expect_binary_val_deriv_eq(exp_var.val(), base_exp_var1.val(),
                             base_exp_var2, test_var.val(),
                             base_test_var1.val(), base_test_var2);
  expect_binary_val_deriv_eq(exp_var.d_, base_exp_var1.d_,
                             base_exp_var2, test_var.d_,
                             base_test_var1.d_, base_test_var2);
}

template <typename FV>
static inline void expect_binary_val_deriv_eq(
    stan::math::fvar<FV> exp_var, int base_exp_var1,
    stan::math::fvar<FV> base_exp_var2, stan::math::fvar<FV> test_var,
    int base_test_var1, stan::math::fvar<FV> base_test_var2) {
  expect_binary_val_deriv_eq(exp_var.val(), base_exp_var1,
                             base_exp_var2.val(), test_var.val(),
                             base_test_var1, base_test_var2.val());
  expect_binary_val_deriv_eq(exp_var.d_, base_exp_var1,
                             base_exp_var2.d_, test_var.d_,
                             base_test_var1, base_test_var2.d_);
}

template <typename FV>
static inline void expect_binary_val_deriv_eq(
    stan::math::fvar<FV> exp_var, stan::math::fvar<FV> base_exp_var1,
    double base_exp_var2, stan::math::fvar<FV> test_var,
    stan::math::fvar<FV> base_test_var1, double base_test_var2) {
  expect_binary_val_deriv_eq(exp_var.val(), base_exp_var1.val(),
                             base_exp_var2, test_var.val(),
                             base_test_var1.val(), base_test_var2);
  expect_binary_val_deriv_eq(exp_var.d_, base_exp_var1.d_,
                             base_exp_var2, test_var.d_,
                             base_test_var1.d_, base_test_var2);
}

template <typename FV>
static inline void expect_binary_val_deriv_eq(
    stan::math::fvar<FV> exp_var, double base_exp_var1,
    stan::math::fvar<FV> base_exp_var2, stan::math::fvar<FV> test_var,
    double base_test_var1, stan::math::fvar<FV> base_test_var2) {
  expect_binary_val_deriv_eq(exp_var.val(), base_exp_var1,
                             base_exp_var2.val(), test_var.val(),
                             base_test_var1, base_test_var2.val());
  expect_binary_val_deriv_eq(exp_var.d_, base_exp_var1,
                             base_exp_var2.d_, test_var.d_,
                             base_test_var1, base_test_var2.d_);
}
#endif
