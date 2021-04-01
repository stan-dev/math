#ifndef STAN_MATH_REV_FUN_AS_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_REV_FUN_AS_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Converts a `var_value<T>` with inner Eigen matrix type to an `var_value<T>`
 *  with an inner array.
 *
 * @tparam T Type of `var_value<T>` with inner `Eigen::Matrix` or expression
 * @param v Specified `var_value<T>` with inner `Eigen::Matrix` or expression.
 * @return `var_value<>` with Matrix converted to an array.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto as_array_or_scalar(T&& v) {
  return v.array();
}

}  // namespace math
}  // namespace stan

#endif
