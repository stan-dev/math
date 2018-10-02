#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_SCALAR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_rev_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_val_deriv_eq.hpp>
#include <vector>

template <typename F>
void expect_rev_scalar_value() {
  using stan::math::var;
  using std::vector;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<var> y = build_rev_vector<F>();
    vector<var> z = build_rev_vector<F>();
    var fz = F::template apply<var>(z[i]);
    expect_val_deriv_eq(F::apply_base(y[i]), y[i], fz, z[i]);
  }
}

#endif
