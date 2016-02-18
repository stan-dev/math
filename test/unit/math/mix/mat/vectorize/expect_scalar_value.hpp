#ifndef TEST_UNIT_MATH_MIX_MAT_TECTORIZE_EXPECT_SCALAR_TALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_TECTORIZE_EXPECT_SCALAR_TALUE_HPP

#include <vector>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <test/unit/math/mix/mat/vectorize/build_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_fvar_var_eq.hpp>

template <typename F, typename T>
void expect_scalar_value() {
  using stan::test::expect_match_return_t;
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<T> y = build_vector<F>(vector<T>(), i);
    vector<T> z = build_vector<F>(vector<T>(), i);
    T fz = F::template apply<T>(z[i]);
    expect_fvar_var_eq(F::apply_base(y[i]), y[i], fz, z[i]);
  }
  expect_match_return_t<F, T, T>();
  expect_match_return_t<F, std::vector<T>, std::vector<T> >();
}

#endif
