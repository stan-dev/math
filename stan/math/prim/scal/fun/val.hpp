#ifndef STAN_MATH_PRIM_SCAL_FUN_VAL_HPP
#define STAN_MATH_PRIM_SCAL_FUN_VAL_HPP

#include <stan/math/prim/scal/meta/is_fr_var.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

// helper descent function arithmetic case
template <class T, std::enable_if_t<
 std::is_arithmetic<T>::value>* =nullptr>
inline auto
val_helper (T const& t) {
 return t;
}

// helper descent function for stan object
template <class T, std::enable_if_t<
 !std::is_arithmetic<T>::value>* =nullptr>
inline auto
val_helper (T const& t) {
 return t.val();
}

/**
 * Returns the underlying value of the
 * stan type, zero or multiple levels deep
 * if required by template parameter S.
 *
 * This is the termination condition
 * of the function template.
 *
 * @tparam T type argument
 * @tparam S return type sought
 * @param t argument
 * @param s optional return type arg
 * @return underlying value
 */
template <class T, class S = decltype(val_helper(T())),
 std::enable_if_t<std::is_same<T, S>::value>* =nullptr>
inline S
val(T const& t, S const& = S()) {
  return t;  // terminate descent regardless of type T
}

/**
 * Returns the underlying value of the
 * (f)var type, zero or multiple levels deep
 * if required by template parameter S.
 *
 * @tparam T type argument
 * @tparam S return type sought
 * @param t argument.
 * @param s optional return type arg
 * @return value of variable.
 */
template <class T, class S = decltype(val_helper(T())),
 std::enable_if_t<!std::is_same<T, S>::value>* =nullptr>
inline S
val(T const& t, S const& s = S()) {
  // fundamental types require
  // val to be fully qualified
  return stan::math::val(val_helper(t), s);
}

/**
 * Returns the underlying value of the
 * stan type, multiple levels deep if
 * required by template parameter S.
 *
 * @tparam T complex type argument
 * @tparam S complex return type sought
 * @param t argument
 * @param s optional return type arg
 * @return underlying value
 */
template <class T, class S = decltype(val_helper(T()))>
inline std::complex<S>
val(std::complex<T> const& t, S const& s = S()) {
  return std::complex<S>(
      // fundamental types require
      // val to be fully qualified
      stan::math::val(t.real(), s),
      stan::math::val(t.imag(), s));
}

}  // namespace math
}  // namespace stan
#endif
