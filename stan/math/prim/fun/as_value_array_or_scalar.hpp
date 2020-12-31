#ifndef STAN_MATH_PRIM_FUN_AS_VALUE_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_FUN_AS_VALUE_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Extract value from scalar types.
 *
 * @tparam T Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto as_value_array_or_scalar(T&& v) {
  return value_of(std::forward<T>(v));
}

/**
 * Extract the values from a Matrix and convert to an array.
 *
 * @tparam T Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename T, require_matrix_t<T>* = nullptr,
          require_not_st_arithmetic<T>* = nullptr>
inline auto as_value_array_or_scalar(T&& v) {
  return make_holder(
      [](auto&& x) { return value_of(std::forward<decltype(x)>(x)); },
      std::forward<T>(v));
}

/**
 * Extract the values from a Matrix and convert to an array.
 *
 * @tparam T Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename T, require_matrix_t<T>* = nullptr,
          require_st_arithmetic<T>* = nullptr>
inline auto as_value_array_or_scalar(T&& v) {
  return as_array_or_scalar(std::forward<T>(v));
}

/**
 * Converts a std::vector type to an array.
 *
 * @tparam T Type of scalar element.
 * @param v Specified vector.
 * @return Matrix converted to an array.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto as_value_array_or_scalar(T&& v) {
  using T_map
      = Eigen::Map<const Eigen::Array<value_type_t<T>, Eigen::Dynamic, 1>>;
  return make_holder(
      [](auto&& x) { return value_of(T_map(x.data(), x.size())); },
      std::forward<T>(v));
}

}  // namespace math
}  // namespace stan

#endif
