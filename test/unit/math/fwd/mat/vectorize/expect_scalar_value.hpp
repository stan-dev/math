#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_SCALAR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <vector>
#include <test/unit/math/fwd/mat/vectorize/build_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_eq.hpp>

template <typename F, typename T>
void expect_scalar_value() {
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<T> y = 
      build_vector<F>(vector<T>(), i);
    T fy = F::template apply<T>(y[i]);
    T exp_y = F::apply_base(y[i]);
    expect_eq(exp_y, fy);
  }
}

#endif
