#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_eq.hpp>

template <typename F>
void expect_scalar_value() {
  using stan::math::var;
  using std::vector;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<var> y = build_vector<F>();
    var fy = F::template apply<var>(y[i]);

    vector<var> z = build_vector<F>();
    var fz = F::apply_base(z[i]);

    expect_eq(fz, z[i], fy, y[i]);
  }
}

#endif
