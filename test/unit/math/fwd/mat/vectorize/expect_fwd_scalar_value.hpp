#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_SCALAR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_SCALAR_VALUE_HPP

#include <stan/math/fwd/mat.hpp>
#include <vector>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F, typename T>
void expect_fwd_scalar_value() {
  using std::vector;
  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<T> y = 
      build_fwd_vector<F>(vector<T>(), i);
    T fy = F::template apply<T>(y[i]);
    T exp_y = F::apply_base(y[i]);
    expect_val_deriv_eq(exp_y, fy);
  }
}

#endif
