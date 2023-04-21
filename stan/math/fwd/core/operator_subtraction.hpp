#ifndef STAN_MATH_FWD_CORE_OPERATOR_SUBTRACTION_HPP
#define STAN_MATH_FWD_CORE_OPERATOR_SUBTRACTION_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/core/std_complex.hpp>
#include <stan/math/prim/core/operator_subtraction.hpp>

namespace stan {
namespace math {

/**
 * Return the difference of the specified arguments.
 *
 * @tparam T type of values and tangents
 * @param x1 first argument
 * @param x2 second argument
 * @return difference of the arguments
 */
template <typename T>
inline fvar<T> operator-(const fvar<T>& x1, const fvar<T>& x2) {
  return fvar<T>(x1.val_ - x2.val_, x1.d_ - x2.d_);
}

/**
 * Return the difference of the specified arguments.
 *
 * @tparam T type of values and tangents
 * @param x1 first argument
 * @param x2 second argument
 * @return difference of the arguments
 */
template <typename T>
inline fvar<T> operator-(double x1, const fvar<T>& x2) {
  return fvar<T>(x1 - x2.val_, -x2.d_);
}

/**
 * Return the difference of the specified arguments.
 *
 * @tparam T type of values and tangents
 * @param x1 first argument
 * @param x2 second argument
 * @return difference of the arguments
 */
template <typename T>
inline fvar<T> operator-(const fvar<T>& x1, double x2) {
  return fvar<T>(x1.val_ - x2, x1.d_);
}

/**
 * Return the difference of the two complex fvar<T> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return subtraction of arguments
 */
template <typename T>
inline std::complex<stan::math::fvar<T>>
operator-(const std::complex<stan::math::fvar<T>>& x,
          const std::complex<stan::math::fvar<T>>& y) {
  return internal::complex_subtract(x, y);
}

/**
 * Return the difference of std::complex<double> and
 * std::complex<fvar<T>> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return subtraction of arguments
 */
template <typename T>
inline std::complex<stan::math::fvar<T>>
operator-(const std::complex<double>& x,
          const std::complex<stan::math::fvar<T>>& y) {
  return internal::complex_subtract(x, y);
}

/**
 * Return the difference of std::complex<double> and
 * std::complex<fvar<T>> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return subtraction of arguments
 */
template <typename T>
inline std::complex<stan::math::fvar<T>>
operator-(const std::complex<stan::math::fvar<T>>& x,
          const std::complex<double>& y) {
  return internal::complex_subtract(x, y);
}

}  // namespace math
}  // namespace stan
#endif
