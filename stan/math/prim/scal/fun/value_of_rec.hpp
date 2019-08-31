#ifndef STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_REC_HPP

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
 * @tparam T Type of scalar.
 * @param x Scalar to convert to double.
 * @return Value of scalar cast to a double.
 */
template <typename T, require_floating_point<T>...>
inline auto&& value_of_rec(T&& x) {
  return std::forward<T>(x);
}

template <typename T, require_arithmetic<T>...,
          require_not_floating_point<T>...>
inline auto value_of_rec(T&& x) {
  return static_cast<double>(x);
}

}  // namespace math
}  // namespace stan
#endif
