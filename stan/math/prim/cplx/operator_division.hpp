#ifndef STAN_MATH_PRIM_CPLX_OPERATOR_DIVISION_HPP
#define STAN_MATH_PRIM_CPLX_OPERATOR_DIVISION_HPP

#include <stan/math/prim/scal/fun/copysign.hpp>
#include <stan/math/prim/scal/fun/complex_promote.hpp>
#include <stan/math/prim/scal/fun/isfinite.hpp>
#include <stan/math/prim/scal/meta/is_arith_like.hpp>
#include <stan/math/prim/scal/meta/is_complex.hpp>
#include <stan/math/prim/scal/meta/is_fr_var.hpp>
#include <complex>
#include <cmath>

namespace stan {
namespace math {
 
/**
 * Return the result of dividing the first argument by the second.
 *
 * @tparam T type of complex first argument
 * @tparam U type of second argument
 * @param t complex first argument
 * @param u second argument
 * @return first argument divided by second argument
 */
template <class T, class U, std::enable_if_t<is_fr_var<T, U>::value
 && (is_complex<U>::value || is_arith_like<U>::value)>* = nullptr>
inline auto operator/(std::complex<T> const& t, U const& u) {
  return complex_promote<T, U>(t) /= u;
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * @tparam T type of first argument
 * @tparam U type of complex second argument
 * @param t first argument
 * @param u complex second argument
 * @return first argument divided by second argument
 */
template <class T, class U, std::enable_if_t<is_fr_var<T, U>::value
 && is_arith_like<T>::value>* = nullptr>
inline auto operator/(T const& t, std::complex<U> const& u) {
  return complex_promote<T, U>(t) /= u;
}

/**
 * Return the result of dividing the first argument by the second.
 * This function exists because libc++ uses logb and scalbn for
 * complex division, which aren't defined for fvars. Four tightly
 * scoped high priority operator/ overloads for fvar and var
 * forward to this function for implementation.
 *
 * @tparam T type of complex value
 * @param t first argument
 * @param u second argument
 * @return first argument divided by second argument
 */
template <class T>
inline auto
operator_division(std::complex<T> const& t,
                  std::complex<T> const& u) {
  using std::pow;
  T const n(pow(u.real(), 2) + pow(u.imag(), 2));
  T const r((t.real() * u.real() + t.imag() * u.imag()) / n);
  T const i((t.imag() * u.real() - t.real() * u.imag()) / n);
  return std::complex<rm_zeroing_t<T>>(r, i);
}

/**
 * Return the result of dividing the first argument by the second.
 *
 * This overload exists because libc++ uses logb and scalbn for
 * complex division, which aren't defined for fvars.
 *
 * @tparam T type of complex fvar value and tangent
 * @param t first argument
 * @param u second argument
 * @return first argument divided by second argument
 */
template <class T>
inline std::complex<rm_zeroing_t<T>> operator/(
    std::complex<zeroing<T>> const& t,
    std::complex<zeroing<T>> const& u) {
  return stan::math::operator_division(t, u);  // no recursion
}

}  // namespace math
}  // namespace stan
#endif
