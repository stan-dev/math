#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_VALUE_HPP

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <gtest/gtest.h>
#include <vector>

template <typename F>
void expect_prim_binary_scalar_value() {
  using std::vector;
  using stan::math::is_nan;
  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();
  vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  //int, int
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    double v = F::template apply<double>(int_valid_inputs1[i],
    int_valid_inputs2[i]);
    double exp_v = F::apply_base(int_valid_inputs1[i], 
    int_valid_inputs2[i]);
    if (is_nan(v) && is_nan(exp_v)) continue;
    EXPECT_FLOAT_EQ(exp_v, v);
  }
  //double, double
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    double v = F::template apply<double>(valid_inputs1[i], 
    valid_inputs2[i]);
    double exp_v = F::apply_base(valid_inputs1[i], valid_inputs2[i]);
    if (is_nan(v) && is_nan(exp_v)) continue;
    EXPECT_FLOAT_EQ(exp_v, v);
  }
  //int, double
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    double v = F::template apply<double>(int_valid_inputs1[i], 
    valid_inputs2[i]);
    double exp_v = F::apply_base(int_valid_inputs1[i], valid_inputs2[i]);
    if (is_nan(v) && is_nan(exp_v)) continue; 
    EXPECT_FLOAT_EQ(exp_v, v);
  }
  //double, int
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    double v = F::template apply<double>(valid_inputs1[i], 
    int_valid_inputs2[i]);
    double exp_v = F::apply_base(valid_inputs1[i], int_valid_inputs2[i]);
    if (is_nan(v) && is_nan(exp_v)) continue; 
    EXPECT_FLOAT_EQ(exp_v, v);
  }
}

#endif
