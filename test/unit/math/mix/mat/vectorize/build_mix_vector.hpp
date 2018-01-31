#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_BUILD_MIX_VECTOR_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_BUILD_MIX_VECTOR_HPP

#include <stan/math/mix/mat.hpp>
#include <vector>

template <typename F>
static inline std::vector<stan::math::var> build_mix_vector(
    std::vector<stan::math::var> fvar_vector, int seed_index) {
  using stan::math::var;
  using std::vector;
  vector<double> inputs = F::valid_inputs();
  for (size_t i = 0; i < inputs.size(); ++i)
    fvar_vector.push_back(inputs[i]);
  return fvar_vector;
}

template <typename F, typename T>
static inline std::vector<stan::math::fvar<T> > build_mix_vector(
    std::vector<stan::math::fvar<T> > fvar_vector, int seed_index = -1) {
  using stan::math::fvar;
  using std::vector;

  vector<T> val_vector = build_mix_vector<F>(vector<T>(), seed_index);
  vector<T> d_vector;
  if (seed_index != -1)
    d_vector = build_mix_vector<F>(vector<T>(), seed_index);

  for (size_t i = 0; i < val_vector.size(); ++i) {
    // For fvar<fvar<var> >, this ensures that when
    // we test the var, we don't autodiff the same var twice
    if (seed_index == static_cast<int>(i))
      fvar_vector.push_back(fvar<T>(val_vector[i], d_vector[i]));
    else
      fvar_vector.push_back(fvar<T>(val_vector[i]));
  }
  return fvar_vector;
}

#endif
