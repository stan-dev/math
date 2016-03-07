#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_SCALAR_VALUE_HPP

#include <vector>
#include <test/unit/math/mix/mat/vectorize/build_mix_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F, typename T>
void expect_mix_scalar_value() {
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<T> y = build_mix_vector<F>(vector<T>(), i);
    vector<T> z = build_mix_vector<F>(vector<T>(), i);
    T fz = F::template apply<T>(z[i]);
    expect_val_deriv_eq(F::apply_base(y[i]), y[i], fz, z[i]);
  }
}

#endif
