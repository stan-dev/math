#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_SCALAR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <vector>

template <typename F>
void expect_prim_scalar_value() {
  using std::vector;
  vector<double> valid_inputs = F::valid_inputs();
  vector<int> int_valid_inputs = F::int_valid_inputs();
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int input = int_valid_inputs[i];
    double exp_v = F::apply_base(static_cast<double>(input));
    double v = F::template apply<double>(input);
    expect_val_eq(exp_v, v);
  }
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    expect_val_eq(F::apply_base(valid_inputs[i]),
                  F::template apply<double>(valid_inputs[i]));
  }
}

#endif
