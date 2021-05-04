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
 * Extract the value from an object. For eigen types and `std::vectors`
 * convert to an eigen array and for scalars return a scalar.
 * @tparam T A stan scalar, eigen vector, or `std::vector`.
 * @param v Specified value.
 * @return Same value.
 */
template <typename T>
inline auto as_value_array_or_scalar(T&& v) {
  return value_of(as_array_or_scalar(std::forward<T>(v)));
}

}  // namespace math
}  // namespace stan

#endif
