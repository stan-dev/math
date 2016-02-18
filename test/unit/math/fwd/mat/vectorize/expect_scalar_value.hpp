#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <vector>
#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fvar_eq.hpp>

template <typename F, typename T>
void expect_scalar_value() {
  using stan::test::expect_match_return_t;
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<T> y = 
      build_vector<F>(vector<T>(), i);
    T fy = F::template apply<T>(y[i]);
    T exp_y = F::apply_base(y[i]);
    expect_fvar_eq(exp_y, fy);
  }
  expect_match_return_t<F, T, T>();
  expect_match_return_t<F, std::vector<T>, std::vector<T> >();

}

#endif
