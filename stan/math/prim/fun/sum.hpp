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
 * @param v Specified value.
 * @return Same value (the sum of one value).
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T sum(T v) {
  return v;
}

/**
 * Return the sum of the values in the specified standard vector.
 *
 * @tparam T Type of elements summed.
 * @param xs Standard vector to sum.
 * @return Sum of elements.
 */
template <typename T>
inline T sum(const std::vector<T>& xs) {
  return std::accumulate(xs.begin(), xs.end(), T{0});
}

/**
 * Returns the sum of the coefficients of the specified
 * Eigen Matrix, Array or expression.
 *
 * @tparam T Type of argument
 * @param v argument
 * @return Sum of coefficients of argument.
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline value_type_t<T> sum(const T& v) {
  return v.sum();
}

}  // namespace math
}  // namespace stan

#endif
