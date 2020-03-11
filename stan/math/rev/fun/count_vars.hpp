#ifndef STAN_MATH_REV_FUN_COUNT_VARS_HPP
#define STAN_MATH_REV_FUN_COUNT_VARS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

inline size_t count_vars_impl(size_t count) { return count; }

template <typename... Pargs>
size_t count_vars_impl(size_t count, const std::vector<var>& x,
                       const Pargs&... args);

template <typename T, require_t<is_var<scalar_type_t<T>>>..., typename... Pargs>
size_t count_vars_impl(size_t count, const std::vector<T>& x,
                       const Pargs&... args);

template <typename Mat, require_eigen_vt<is_var, Mat>..., typename... Pargs>
size_t count_vars_impl(size_t count, const Mat& x, const Pargs&... args);

template <typename... Pargs>
size_t count_vars_impl(size_t count, const var& x, const Pargs&... args);

template <typename... Pargs, typename Arith,
          require_arithmetic_t<scalar_type_t<Arith>>...>
size_t count_vars_impl(size_t count, Arith& x, const Pargs&... args);

template <typename... Pargs>
size_t count_vars_impl(size_t count, const std::vector<var>& x,
                       const Pargs&... args) {
  return count_vars_impl(count + x.size(), args...);
}

template <typename T, require_t<is_var<scalar_type_t<T>>>..., typename... Pargs>
size_t count_vars_impl(size_t count, const std::vector<T>& x,
                       const Pargs&... args) {
  for (size_t i = 0; i < x.size(); i++) {
    count = count_vars_impl(count, x[i]);
  }
  return count_vars_impl(count, args...);
}

template <typename Mat, require_eigen_vt<is_var, Mat>..., typename... Pargs>
size_t count_vars_impl(size_t count, const Mat& x, const Pargs&... args) {
  return count_vars_impl(count + x.size(), args...);
}

template <typename... Pargs>
size_t count_vars_impl(size_t count, const var& x, const Pargs&... args) {
  return count_vars_impl(count + 1, args...);
}

template <typename... Pargs, typename Arith,
          require_arithmetic_t<scalar_type_t<Arith>>...>
size_t count_vars_impl(size_t count, Arith& x, const Pargs&... args) {
  return count_vars_impl(count, args...);
}

/**
 * Count the number of scalars of type T in the input argument list
 *
 * @tparam Pargs Types of input arguments
 * @return Number of scalars of type T in input
 */
template <typename... Pargs>
size_t count_vars(const Pargs&... args) {
  return count_vars_impl(0, args...);
}

}  // namespace internal

}  // namespace math
}  // namespace stan
#endif
