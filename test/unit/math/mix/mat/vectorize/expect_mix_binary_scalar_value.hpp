#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_SCALAR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>

template <typename F, typename FV>
void expect_mix_binary_scalar_value() {
  using stan::math::fvar;
  using std::vector;

  vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();
  // FV, int
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    int input1 = int_valid_inputs1[i];
    int input2 = int_valid_inputs2[i];
    FV y1 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV y2 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV z1 = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV z2 = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV fz = F::template apply<FV>(y2, input2);
    expect_binary_val_deriv_eq(F::apply_base(y1, input2), y1, input2, fz, y2,
                               input2);
    fz = F::template apply<FV>(input1, z2);
    expect_binary_val_deriv_eq(F::apply_base(input1, z1), input1, z1, fz,
                               input1, z2);
  }
  // FV, double
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    double input1 = valid_inputs1[i];
    double input2 = valid_inputs2[i];
    FV y1 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV y2 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV z1 = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV z2 = build_binary_vector2<F>(vector<FV>(), i)[i];
    FV fz = F::template apply<FV>(y2, input2);
    expect_binary_val_deriv_eq(F::apply_base(y1, input2), y1, input2, fz, y2,
                               input2);
    fz = F::template apply<FV>(input1, z2);
    expect_binary_val_deriv_eq(F::apply_base(input1, z1), input1, z1, fz,
                               input1, z2);
  }
  // FV, FV
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    FV y1 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV y2 = build_binary_vector1<F>(vector<FV>(), i)[i];
    FV z1 = build_binary_vector2<F>(vector<FV>())[i];
    FV z2 = build_binary_vector2<F>(vector<FV>())[i];
    FV fz = F::template apply<FV>(y2, z2);
    expect_binary_val_deriv_eq(F::apply_base(y1, z1), y1, z1, fz, y2, z2);
    y1 = build_binary_vector1<F>(vector<FV>())[i];
    y2 = build_binary_vector1<F>(vector<FV>())[i];
    z1 = build_binary_vector2<F>(vector<FV>(), i)[i];
    z2 = build_binary_vector2<F>(vector<FV>(), i)[i];
    fz = F::template apply<FV>(y2, z2);
    expect_binary_val_deriv_eq(F::apply_base(y1, z1), y1, z1, fz, y2, z2);
    y1 = build_binary_vector1<F>(vector<FV>(), i)[i];
    y2 = build_binary_vector1<F>(vector<FV>(), i)[i];
    z1 = build_binary_vector2<F>(vector<FV>(), i)[i];
    z2 = build_binary_vector2<F>(vector<FV>(), i)[i];
    fz = F::template apply<FV>(y2, z2);
    expect_binary_val_deriv_eq(F::apply_base(y1, z1), y1, z1, fz, y2, z2);
  }
}

#endif
