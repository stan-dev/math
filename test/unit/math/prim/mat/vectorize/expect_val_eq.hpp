#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_EXPECT_VAL_EQ_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_EXPECT_VAL_EQ_HPP

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <gtest/gtest.h>

void expect_val_eq(double exp_val, double test_val) {
  using stan::math::is_nan;
  if (is_nan(exp_val) && is_nan(test_val)) {
    SUCCEED();
    return;
  }
  EXPECT_FLOAT_EQ(exp_val, test_val);
}

#endif
