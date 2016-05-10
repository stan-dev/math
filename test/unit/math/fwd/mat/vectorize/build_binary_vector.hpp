#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_BINARY_VECTOR_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_BINARY_VECTOR_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/prim/mat/vectorize/build_binary_vector.hpp>
#include <vector>

template <typename F, typename FV>
static inline std::vector<stan::math::fvar<FV> >
build_binary_vector1(std::vector<stan::math::fvar<FV> > fvar_vector,
int seed_index = -1) { 
  using std::vector;
  using stan::math::fvar;

  vector<FV> template_v = build_binary_vector1<F>(vector<FV>(), seed_index);
  vector<fvar<FV> > result_v;

  for (size_t i = 0; i < template_v.size(); ++i) {
  
    // For fvar<fvar<double> >, this will fill in 
    // all four components
    if (seed_index == static_cast<int>(i))
      result_v.push_back(fvar<FV>(template_v[i], template_v[i]));
    else
      result_v.push_back(fvar<FV>(template_v[i]));
  }
  return result_v;
}

template <typename F, typename FV>
static inline std::vector<stan::math::fvar<FV> >
build_binary_vector2(std::vector<stan::math::fvar<FV> > fvar_vector,
int seed_index = -1) { 
  using std::vector;
  using stan::math::fvar;

  vector<FV> template_v = build_binary_vector2<F>(vector<FV>(), seed_index);
  vector<fvar<FV> > result_v;
  for (size_t i = 0; i < template_v.size(); ++i) {
  
    // For fvar<fvar<double> >, this will fill in 
    // all four components
    if (seed_index == static_cast<int>(i))
      result_v.push_back(fvar<FV>(template_v[i], template_v[i]));
    else
      result_v.push_back(fvar<FV>(template_v[i]));
  }
  return result_v;
}

#endif
