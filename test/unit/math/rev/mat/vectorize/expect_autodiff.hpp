#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_AUTODIFF_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_AUTODIFF_HPP

#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>

template <typename F>
static inline void expect_autodiff(double test_val, double test_adj) {
  using stan::math::var;

  var x = test_val;
  F::apply_base(x).grad();
  EXPECT_FLOAT_EQ(x.adj(), test_adj);
}

#endif
