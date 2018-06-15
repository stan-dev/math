#ifndef STAN_MATH_PRIM_SCAL_META_IS_ARITH_LIKE_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_ARITH_LIKE_HPP

#include <stan/math/prim/scal/meta/is_fr_var.hpp>
#include <stan/math/prim/scal/meta/is_complex.hpp>
#include <type_traits>

namespace stan {

/**
 * std::true_type if T is fvar, var, or arithmetic
 * std::false_type if T is complex or other
 * restriction on complex matches std::is_arithmetic
 */
template <class T>
struct is_arith_like_helper
 : std::integral_constant<bool,
    !is_complex<T>::value
    && (is_fr_var<T>::value
     || std::is_arithmetic<T>::value)> {};

/**
 * std::true_type if any parameter is fvar, var, 
 * or arithmetic. Restriction on complex matches
 * std::is_arithmetic.
 */
template <class T, class U = double, class V = double,
 class W = double, class X = double, class Y = double>
struct is_arith_like :
 std::integral_constant<bool,
  is_arith_like_helper<T>::value
  || is_arith_like_helper<U>::value
  || is_arith_like_helper<V>::value
  || is_arith_like_helper<W>::value
  || is_arith_like_helper<X>::value
  || is_arith_like_helper<Y>::value> {};

}  // namespace stan
#endif
