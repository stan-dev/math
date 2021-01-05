#ifndef STAN_MATH_REV_FUN_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_REV_FUN_AS_COLUMN_VECTOR_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>

namespace stan {
namespace math {

/**
 * Converts input argument to a `var_value<>` with column vector or a scalar.
 * For column vector inputs this is an identity function.
 *
 * @tparam T Type of `var_value`.
 * @param a Specified vector.
 * @return Same vector.
 */
template <typename T, require_var_col_vector_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  return std::forward<T>(a);
}

/**
 * Converts `var_value` with row vector inner type to a `var_value<>`
 * with inner column vector type
 * @tparam T A `var_value<>` with an inner row vector type.
 * @param a Specified vector.
 * @return Transposed vector.
 */
template <typename T, require_var_row_vector_t<T>* = nullptr,
          require_not_var_col_vector_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  return a.transpose();
}

}  // namespace math
}  // namespace stan

#endif
