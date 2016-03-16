#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_SCALAR_VALUE_HPP

#include <gtest/gtest.h>
#include <vector>

template <typename F>
void expect_prim_scalar_value() {
  using std::vector;
  vector<double> valid_inputs = F::valid_inputs();
  vector<int> int_valid_inputs = F::int_valid_inputs();
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int input = int_valid_inputs[i];
    double v = F::template apply<double>(input);
    EXPECT_FLOAT_EQ(F::apply_base(double(input)), v);
  }
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    EXPECT_FLOAT_EQ(F::apply_base(valid_inputs[i]),
                    F::template apply<double>(valid_inputs[i]));
  }
}

#endif
