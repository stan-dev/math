#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_autodiff.hpp>

template <typename F>
void expect_scalar_value() {
  using stan::math::var;
  using std::vector;

  vector<var> y = build_vector<F>();
  for (size_t i = 0; i < y.size(); ++i) {
    var fy = F::template apply<var>(y[i]);

    EXPECT_FLOAT_EQ(F::apply_base(y[i]).val(), fy.val());

    fy.grad();
    expect_autodiff<F>(y[i].val(), y[i].adj());
  }
}

#endif
