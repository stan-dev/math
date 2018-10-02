#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_FWD_VECTOR_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_FWD_VECTOR_HPP

#include <stan/math/fwd/mat.hpp>
#include <vector>

template <typename F>
static inline std::vector<double> build_fwd_vector(
    std::vector<double> double_vector, int seed_index = -1) {
  return F::valid_inputs();
}

template <typename F, typename T>
static inline std::vector<stan::math::fvar<T> > build_fwd_vector(
    std::vector<stan::math::fvar<T> > fvar_vector, int seed_index = -1) {
  using stan::math::fvar;
  using std::vector;

  vector<T> template_v = build_fwd_vector<F>(vector<T>(), seed_index);

  for (size_t i = 0; i < template_v.size(); ++i) {
    // For fvar<fvar<double> >, this will fill in
    // all four components
    if (seed_index == static_cast<int>(i))
      fvar_vector.push_back(fvar<T>(template_v[i], template_v[i]));
    else
      fvar_vector.push_back(fvar<T>(template_v[i]));
  }
  return fvar_vector;
}

#endif
