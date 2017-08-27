#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_BINARY_VECTOR_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_BINARY_VECTOR_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/prim/mat/vectorize/build_binary_vector.hpp>
#include <vector>

template <typename F> static inline
std::vector<stan::math::var> build_binary_vector1(
    std::vector<stan::math::var> template_v, int seed = -1) {
  using std::vector;
  using stan::math::var;
  vector<double> inputs = F::valid_inputs1();
  return vector<var>(inputs.begin(), inputs.end());
}

template <typename F> static inline
std::vector<stan::math::var> build_binary_vector2(
    std::vector<stan::math::var> template_v, int seed = -1) {
  using std::vector;
  using stan::math::var;
  vector<double> inputs = F::valid_inputs2();
  return vector<var>(inputs.begin(), inputs.end());
}

#endif
