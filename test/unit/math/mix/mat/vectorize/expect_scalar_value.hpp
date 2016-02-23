#ifndef TEST_UNIT_MATH_MIX_MAT_TECTORIZE_EXPECT_SCALAR_TALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_TECTORIZE_EXPECT_SCALAR_TALUE_HPP

#include <vector>
#include <test/unit/math/mix/mat/vectorize/build_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_eq.hpp>

template <typename F, typename T>
void expect_scalar_value() {
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<T> y = build_vector<F>(vector<T>(), i);
    vector<T> z = build_vector<F>(vector<T>(), i);
    T fz = F::template apply<T>(z[i]);
    expect_eq(F::apply_base(y[i]), y[i], fz, z[i]);
  }
}

#endif
