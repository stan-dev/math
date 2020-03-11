#ifndef STAN_MATH_REV_FUN_ACCUMULATE_ADJOINTS_HPP
#define STAN_MATH_REV_FUN_ACCUMULATE_ADJOINTS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

inline double* accumulate_adjoints(double* x) { return x; }

template <typename... Pargs>
double* accumulate_adjoints(double* dest, const var& x, const Pargs&... args);

template <typename... Pargs>
double* accumulate_adjoints(double* dest, const std::vector<var>& x,
                            const Pargs&... args);

template <typename T, require_t<is_var<scalar_type_t<T>>>..., typename... Pargs>
double* accumulate_adjoints(double* dest, const std::vector<T>& x,
                            const Pargs&... args);

template <typename... Pargs, typename Mat, require_eigen_vt<is_var, Mat>...>
double* accumulate_adjoints(double* dest, const Mat& x, const Pargs&... args);

template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>>...,
          typename... Pargs>
double* accumulate_adjoints(double* dest, Arith&& x, const Pargs&... args);

template <typename... Pargs>
double* accumulate_adjoints(double* dest, const var& x, const Pargs&... args) {
  *dest += x.adj();
  return accumulate_adjoints(dest + 1, args...);
}

template <typename... Pargs>
double* accumulate_adjoints(double* dest, const std::vector<var>& x,
                            const Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest[i] += x[i].adj();
  }
  return accumulate_adjoints(dest + x.size(), args...);
}

template <typename T, require_t<is_var<scalar_type_t<T>>>..., typename... Pargs>
double* accumulate_adjoints(double* dest, const std::vector<T>& x,
                            const Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest = accumulate_adjoints(dest, x[i]);
  }
  return accumulate_adjoints(dest, args...);
}

template <typename... Pargs, typename Mat, require_eigen_vt<is_var, Mat>...>
double* accumulate_adjoints(double* dest, const Mat& x, const Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest[i] += x(i).adj();
  }
  return accumulate_adjoints(dest + x.size(), args...);
}

template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>>...,
          typename... Pargs>
double* accumulate_adjoints(double* dest, Arith&& x, const Pargs&... args) {
  return accumulate_adjoints(dest, args...);
}

}  // namespace internal

}  // namespace math
}  // namespace stan
#endif
