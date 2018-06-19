#ifndef STAN_MATH_PRIM_SCAL_META_IS_FR_VAR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_FR_VAR_HPP

#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/rm_complex.hpp>
#include <type_traits>

namespace stan {

/**
 * Is std::true_type if T is fvar or var.
 * Strips outer complex types.
 * Note rm_complex also does rm_zeroing.
 */
template <class T>
struct is_fr_var_helper
    : std::integral_constant<bool, is_fvar<rm_complex_t<T>>::value
                                       || is_var<rm_complex_t<T>>::value> {};

/**
 * Is std::true_type if any parameter is fvar or var.
 * Strips outer complex types..
 */
template <class T, class U = double, class V = double, class W = double,
          class X = double, class Y = double>
struct is_fr_var
    : std::integral_constant<
          bool, is_fr_var_helper<T>::value || is_fr_var_helper<U>::value
                    || is_fr_var_helper<V>::value || is_fr_var_helper<W>::value
                    || is_fr_var_helper<X>::value
                    || is_fr_var_helper<Y>::value> {};

}  // namespace stan
#endif
