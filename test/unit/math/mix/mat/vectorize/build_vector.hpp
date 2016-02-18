#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_BUILD_VECTOR_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_BUILD_VECTOR_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/rev/core/var.hpp>
#include <vector>

template <typename F>
static inline std::vector<stan::math::var>
build_vector(std::vector<stan::math::var> fvar_vector, int seed_index) {
  using std::vector;
  using stan::math::var;
  vector<double> inputs = F::valid_inputs();
  for (size_t i = 0; i < inputs.size(); ++i)
    fvar_vector.push_back(inputs[i]);
  return fvar_vector;
}

template <typename F, typename T>
static inline std::vector<stan::math::fvar<T> >
build_vector(std::vector<stan::math::fvar<T> > fvar_vector,
                          int seed_index = -1) { 
  using std::vector;
  using stan::math::fvar;

  vector<T> val_vector =
    build_vector<F>(vector<T>(), seed_index);
  vector<T> d_vector; 
  if (seed_index != -1)
    d_vector = build_vector<F>(vector<T>(), seed_index);

  for (size_t i = 0; i < val_vector.size(); ++i) {
    if (seed_index == static_cast<int>(i))
      fvar_vector.push_back(fvar<T>(val_vector[i], d_vector[i]));
    else
      fvar_vector.push_back(fvar<T>(val_vector[i]));
  }
  return fvar_vector;
}

#endif
