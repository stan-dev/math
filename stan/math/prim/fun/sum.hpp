#ifndef STAN_MATH_PRIM_FUN_SUM_HPP
#define STAN_MATH_PRIM_FUN_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <numeric>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns specified input value.
 *
 * @tparam T Type of element.
 * @param m Specified value.
 * @return Same value (the sum of one value).
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T sum(T&& m) {
  return std::forward<T>(m);
}

/**
 * Return the sum of the values in the specified standard vector.
 *
 * @tparam T Type of elements summed.
 * @param m Standard vector to sum.
 * @return Sum of elements.
 */
template <typename T, require_not_var_t<T>* = nullptr>
inline T sum(const std::vector<T>& m) {
  return std::accumulate(m.begin(), m.end(), T{0});
}

/**
 * Returns the sum of the coefficients of the specified
 * Eigen Matrix, Array or expression.
 *
 * @tparam T Type of argument
 * @param m argument
 * @return Sum of coefficients of argument.
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline value_type_t<T> sum(const T& m) {
  return m.sum();
}

/**
 * Returns the sum of the coefficients of the specified
 * Eigen Matrix, Array or expression of complex type.
 *
 * @tparam T Type of argument
 * @param m argument
 * @return Sum of coefficients of argument.
 */
template <typename T, require_eigen_vt<is_complex, T>* = nullptr>
inline value_type_t<T> sum(const T& m) {
  return m.sum();
}

}  // namespace math
}  // namespace stan

#endif
