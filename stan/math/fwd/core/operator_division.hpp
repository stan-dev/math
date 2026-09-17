#ifndef STAN_MATH_FWD_CORE_OPERATOR_DIVISION_HPP
#define STAN_MATH_FWD_CORE_OPERATOR_DIVISION_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/core/operator_division.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the result of dividing the first argument by the second.
 *
 * @tparam T type of fvar value and tangent
 * @param x1 first argument
 * @param x2 second argument
 * @return first argument divided by second argument
 */
template <typename T>
inline fvar<T> operator/(const fvar<T>& x1, const fvar<T>& x2) {
  return fvar<T>(x1.val_ / x2.val_,
                 (x1.d_ * x2.val_ - x1.val_ * x2.d_) / (x2.val_ * x2.val_));
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @tparam T type of fvar value and tangent
 * @param x1 first argument
 * @param x2 second argument
 * @return first argument divided by second argument
 */
template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline fvar<T> operator/(const fvar<T>& x1, U x2) {
  return fvar<T>(x1.val_ / static_cast<double>(x2),
                 x1.d_ / static_cast<double>(x2));
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @tparam T type of fvar value and tangent
 * @param x1 first argument
 * @param x2 second argument
 * @return first argument divided by second argument
 */
template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline fvar<T> operator/(U x1, const fvar<T>& x2) {
  return fvar<T>(static_cast<double>(x1) / x2.val_,
                 -static_cast<double>(x1) * x2.d_ / (x2.val_ * x2.val_));
}

template <typename T>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>>& x1,
                                       const std::complex<fvar<T>>& x2) {
  return internal::complex_divide(x1, x2);
}
template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>>& x1,
                                       const std::complex<U>& x2) {
  return internal::complex_divide(x1, x2);
}
template <typename T>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>>& x1,
                                       const fvar<T>& x2) {
  return internal::complex_divide(x1, x2);
}
template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<fvar<T>>& x1, U x2) {
  return internal::complex_divide(x1, x2);
}

template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<U>& x1,
                                       const std::complex<fvar<T>>& x2) {
  return internal::complex_divide(x1, x2);
}
template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline std::complex<fvar<T>> operator/(const std::complex<U>& x1,
                                       const fvar<T>& x2) {
  return internal::complex_divide(x1, x2);
}

template <typename T>
inline std::complex<fvar<T>> operator/(const fvar<T>& x1,
                                       const std::complex<fvar<T>>& x2) {
  return internal::complex_divide(x1, x2);
}
template <typename T, typename U,
          typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline std::complex<fvar<T>> operator/(const fvar<T>& x1,
                                       const std::complex<U>& x2) {
  return internal::complex_divide(x1, x2);
}

template <typename T, typename U, require_arithmetic_t<U>* = nullptr>
inline std::complex<fvar<T>> operator/(U x1, const std::complex<fvar<T>>& x2) {
  return internal::complex_divide(x1, x2);
}

}  // namespace math
}  // namespace stan
#endif
