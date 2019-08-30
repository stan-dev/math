#ifndef STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the value of the specified scalar argument
 * converted to a double value.
 *
 * <p>See the <code>primitive_value</code> function to
 * extract values without casting to <code>double</code>.
 *
 * <p>This function is meant to cover the primitive types. For
 * types requiring pass-by-reference, this template function
 * should be specialized.
 *
 * @tparam T type of scalar.
 * @param x scalar to convert to double
 * @return value of scalar cast to double
 */
template <typename T, enable_if_arithmetic<std::decay_t<T>>* = nullptr>
inline auto&& value_of(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
