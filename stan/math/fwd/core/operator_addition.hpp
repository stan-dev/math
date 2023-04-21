#ifndef STAN_MATH_FWD_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_FWD_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/core/std_complex.hpp>
#include <stan/math/prim/core/operator_addition.hpp>

namespace stan {
namespace math {

/**
 * Return the sum of the specified forward mode addends.
 *
 * @tparam T type of values and tangents
 * @param x1 first addend
 * @param x2 second addend
 * @return sum of addends
 */
template <typename T>
inline fvar<T> operator+(const fvar<T>& x1, const fvar<T>& x2) {
  return fvar<T>(x1.val_ + x2.val_, x1.d_ + x2.d_);
}

/**
 * Return the sum of the specified double and forward mode addends.
 *
 * @tparam T type of values and tangents
 * @param x1 first addend
 * @param x2 second addend
 * @return sum of addends
 */
template <typename T>
inline fvar<T> operator+(double x1, const fvar<T>& x2) {
  return fvar<T>(x1 + x2.val_, x2.d_);
}

/**
 * Return the sum of the specified forward mode and double addends.
 *
 * @tparam T type of values and tangents
 * @param x1 first addend
 * @param x2 second addend
 * @return sum of addends
 */
template <typename T>
inline fvar<T> operator+(const fvar<T>& x1, double x2) {
  return fvar<T>(x1.val_ + x2, x1.d_);
}

/**
 * Return the sum of the two complex fvar<T> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return sum of arguments
 */
template <typename T>
inline std::complex<stan::math::fvar<T>> operator+(
    const std::complex<stan::math::fvar<T>>& x,
    const std::complex<stan::math::fvar<T>>& y) {
  return internal::complex_add(x, y);
}

/**
 * Return the sum of std::complex<double> and
 * std::complex<fvar<T>> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return sum of arguments
 */
template <typename T>
inline std::complex<stan::math::fvar<T>> operator+(
    const std::complex<double>& x, const std::complex<stan::math::fvar<T>>& y) {
  return internal::complex_add(x, y);
}

/**
 * Return the sum of std::complex<double> and
 * std::complex<fvar<T>> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return sum of arguments
 */
template <typename T>
inline std::complex<stan::math::fvar<T>> operator+(
    const std::complex<stan::math::fvar<T>>& x, const std::complex<double>& y) {
  return internal::complex_add(x, y);
}

}  // namespace math
}  // namespace stan
#endif
