#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_FWD_BINARY_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_FWD_BINARY_SCALAR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <vector>
#include <gtest/gtest.h>

template <typename F, typename FV>
void expect_fwd_binary_scalar_value() {
  using stan::math::fvar;
  using std::vector;

  vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();
  //FV, int
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    int input1 = int_valid_inputs1[i];
    int input2 = int_valid_inputs2[i];
    FV y = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV z = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV fz = F::template apply<FV>(y, input2);
    expect_val_deriv_eq(F::apply_base(y, input2), fz);
    fz = F::template apply<FV>(input1, z);
    expect_val_deriv_eq(F::apply_base(input1, z), fz);
  }
  //FV, double
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    double input1 = valid_inputs1[i];
    double input2 = valid_inputs2[i];
    FV y = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV z = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV fz = F::template apply<FV>(y, input2);
    expect_val_deriv_eq(F::apply_base(y, input2), fz);
    fz = F::template apply<FV>(input1, z);
    expect_val_deriv_eq(F::apply_base(input1, z), fz);
  }
  //FV, FV
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    FV y1 = build_binary_vector1<F>(vector<FV>())[i];
    FV y2 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV z1 = build_binary_vector2<F>(vector<FV>())[i];
    FV z2 = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV fz = F::template apply<FV>(y2, z1);
    expect_val_deriv_eq(F::apply_base(y2, z1), fz);
    fz = F::template apply<FV>(y1, z2);
    expect_val_deriv_eq(F::apply_base(y1, z2), fz);
    fz = F::template apply<FV>(y2, z2);
    expect_val_deriv_eq(F::apply_base(y2, z2), fz);
  }
}

#endif
