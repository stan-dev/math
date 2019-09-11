#ifndef STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/meta.hpp>
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
template <typename T, typename = void>
inline double value_of_rec(const T x) {
  return static_cast<double>(x);
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x value
 * @return input value
 */
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value>...>
inline auto&& value_of_rec(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
