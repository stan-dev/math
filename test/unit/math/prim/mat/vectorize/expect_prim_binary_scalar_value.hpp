#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_VALUE_HPP

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename F>
void expect_prim_binary_scalar_value() {
  using std::vector;
  using stan::math::is_nan;
  vector<double> valid_inputs = F::valid_inputs();
  vector<int> int_valid_inputs = F::int_valid_inputs();
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int input1 = int_valid_inputs[i];
    for (size_t j = 0; i < int_valid_inputs.size(); ++i) {
      int input2 = int_valid_inputs[j];
      double v = F::template apply<double>(input1, input2);
      double exp_v = F::apply_base(input1, input2);
      if (is_nan(v) && is_nan(exp_v)) continue;
      EXPECT_FLOAT_EQ(exp_v, v);
    }
  }
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    double input1 = valid_inputs[i];                      
    for (size_t j = 0; i < valid_inputs.size(); ++i) {
      double input2 = valid_inputs[j];
      double v = F::template apply<double>(input1, input2);
      double exp_v = F::apply_base(input1, input2);
      if (is_nan(v) && is_nan(exp_v)) continue;
      EXPECT_FLOAT_EQ(exp_v, v);
    }
  }
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int int_input = int_valid_inputs[i];
    for (size_t j = 0; i < valid_inputs.size(); ++i) {
      double input = valid_inputs[j];
      double v = F::template apply<double>(int_input, input);
      double exp_v = F::apply_base(int_input, input);
      if (is_nan(v) && is_nan(exp_v)) 
        EXPECT_FLOAT_EQ(exp_v, v);
      v = F::template apply<double>(input, int_input);
      exp_v = F::apply_base(input, int_input);
      if (is_nan(v) && is_nan(exp_v)) 
        EXPECT_FLOAT_EQ(exp_v, v);
    }
  }
}

#endif
