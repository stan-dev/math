#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_REV_VECTOR_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_REV_VECTOR_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>

template <typename F>
static inline std::vector<stan::math::var> build_rev_vector() {
  using stan::math::var;
  using std::vector;
  vector<double> inputs = F::valid_inputs();
  return vector<var>(inputs.begin(), inputs.end());
}

#endif
