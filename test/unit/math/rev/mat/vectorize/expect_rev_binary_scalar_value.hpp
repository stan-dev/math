#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_SCALAR_VALUE_HPP

#include <gtest/gtest.h>
#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>

template <typename F>
void expect_rev_binary_scalar_value() {
  using stan::math::var;
  using std::vector;

  vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();
  // var, double
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    int input1 = int_valid_inputs1[i];
    int input2 = int_valid_inputs2[i];
    var y = valid_inputs1[i];
    var z = valid_inputs2[i];
    var fz = F::template apply<var>(y, input2);
    expect_binary_val_deriv_eq(F::apply_base(y, input2), y, input2,
                               fz, y, input2);
    fz = F::template apply<var>(input1, z);
    expect_binary_val_deriv_eq(F::apply_base(input1, z), input1, z,
                               fz, input1, z);
  }
  // var, double
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    double input1 = valid_inputs1[i];
    double input2 = valid_inputs2[i];
    var y = valid_inputs1[i];
    var z = valid_inputs2[i];
    var fz = F::template apply<var>(y, input2);
    expect_binary_val_deriv_eq(F::apply_base(y, input2), y, input2,
                               fz, y, input2);
    fz = F::template apply<var>(input1, z);
    expect_binary_val_deriv_eq(F::apply_base(input1, z), input1, z,
                               fz, input1, z);
  }
  // var, var
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    var y1 = valid_inputs1[i];
    var y2 = valid_inputs1[i];
    var z1 = valid_inputs2[i];
    var z2 = valid_inputs2[i];
    var fz = F::template apply<var>(y2, z2);
    expect_binary_val_deriv_eq(F::apply_base(y1, z1), y1, z1, fz, y2, z2);
  }
}

#endif
