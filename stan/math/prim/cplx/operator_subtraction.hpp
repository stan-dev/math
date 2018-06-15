#ifndef STAN_MATH_PRIM_CPLX_OPERATOR_SUBTRACTION_HPP
#define STAN_MATH_PRIM_CPLX_OPERATOR_SUBTRACTION_HPP

#include <stan/math/prim/scal/fun/copysign.hpp>
#include <stan/math/prim/scal/fun/complex_promote.hpp>
#include <stan/math/prim/scal/fun/isfinite.hpp>
#include <stan/math/prim/scal/meta/is_arith_like.hpp>
#include <stan/math/prim/scal/meta/is_complex.hpp>
#include <stan/math/prim/scal/meta/is_fr_var.hpp>
#include <complex>

namespace stan {
namespace math {
 
/**
 * Return the difference of the specified arguments
 *
 * @tparam T type of complex first argument
 * @tparam U type of second argument
 * @param t complex first argument
 * @param u second argument
 * @return difference
 */
template <class T, class U, std::enable_if_t<is_fr_var<T, U>::value
 && (is_complex<U>::value || is_arith_like<U>::value)>* = nullptr>
inline auto operator-(std::complex<T> const& t, U const& u) {
  return complex_promote<T, U>(t) -= u;
}

/**
 * Return the difference of the specified arguments
 *
 * @tparam T type of first argument
 * @tparam U type of complex second argument
 * @param t first argument
 * @param u complex second argument
 * @return difference
 */
template <class T, class U, std::enable_if_t<is_fr_var<T, U>::value
 && is_arith_like<T>::value>* = nullptr>
inline auto operator-(T const& t, std::complex<U> const& u) {
  return complex_promote<T, U>(t) -= u;
}

}  // namespace math
}  // namespace stan
#endif
