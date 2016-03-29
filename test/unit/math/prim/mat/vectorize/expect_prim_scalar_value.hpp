#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_SCALAR_VALUE_HPP

#include <gtest/gtest.h>
#include <vector>
#include <stan/math/prim/scal/fun/is_nan.hpp>

template <typename F>
void expect_prim_scalar_value() {
  using std::vector;
  using stan::math::is_nan;
  vector<double> valid_inputs = F::valid_inputs();
  vector<int> int_valid_inputs = F::int_valid_inputs();
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int input = int_valid_inputs[i];
    double exp_v = F::apply_base(double(input));
    double v = F::template apply<double>(input);
    if (is_nan(exp_v) && is_nan(v)) continue;
    EXPECT_FLOAT_EQ(exp_v, v);
  }
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    if (is_nan(F::apply_base(valid_inputs[i]))
        && is_nan(F::template apply<double>(valid_inputs[i])))
    continue;
    EXPECT_FLOAT_EQ(F::apply_base(valid_inputs[i]),
                    F::template apply<double>(valid_inputs[i]));
  }
}

#endif
