#ifndef STAN_MATH_FWD_CORE_OPERATOR_MULTIPLICATION_HPP
#define STAN_MATH_FWD_CORE_OPERATOR_MULTIPLICATION_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/core/operator_multiplication.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the two arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return product of arguments
 */
template <typename T>
inline fvar<T> operator*(const fvar<T>& x, const fvar<T>& y) {
  return fvar<T>(x.val_ * y.val_, x.d_ * y.val_ + x.val_ * y.d_);
}

/**
 * Return the product of the two arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return product of arguments
 */
template <typename T>
inline fvar<T> operator*(double x, const fvar<T>& y) {
  return fvar<T>(x * y.val_, x * y.d_);
}

/**
 * Return the product of the two arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return product of arguments
 */
template <typename T>
inline fvar<T> operator*(const fvar<T>& x, double y) {
  return fvar<T>(x.val_ * y, x.d_ * y);
}

/**
 * Return the product of the two complex fvar<T> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return product of arguments
 */
template <typename T>
inline stan::math::complex<stan::math::fvar<T>> operator*(
    const stan::math::complex<stan::math::fvar<T>>& x,
    const stan::math::complex<stan::math::fvar<T>>& y) {
  return internal::complex_multiply(x, y);
}

/**
 * Return the product of stan::math::complex<double> and
 * stan::math::complex<fvar<T>> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return product of arguments
 */
template <typename T>
inline stan::math::complex<stan::math::fvar<T>> operator*(
    const stan::math::complex<double>& x,
    const stan::math::complex<stan::math::fvar<T>>& y) {
  return internal::complex_multiply(x, y);
}

/**
 * Return the product of stan::math::complex<double> and
 * stan::math::complex<fvar<T>> arguments.
 *
 * @tparam value and tangent type for variables
 * @param[in] x first argument
 * @param[in] y second argument
 * @return product of arguments
 */
template <typename T>
inline stan::math::complex<stan::math::fvar<T>> operator*(
    const stan::math::complex<stan::math::fvar<T>>& x,
    const stan::math::complex<double>& y) {
  return internal::complex_multiply(x, y);
}

}  // namespace math
}  // namespace stan
#endif
