#ifndef STAN_MATH_PRIM_META_INCLUDE_SUMMAND_HPP
#define STAN_MATH_PRIM_META_INCLUDE_SUMMAND_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_constant.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <type_traits>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Template metaprogram to calculate whether a summand
 * needs to be included in a proportional (log) probability
 * calculation.  For usage, the first boolean parameter
 * should be set to <code>true</code> if calculating
 * a term up to proportionality.  Other type parameters
 * should be included for all of the types of variables
 * in a term.
 *
 * The metaprogram can take an arbitrary number of types.
 *
 * The <code>value</code> enum will be <code>true</code> if the
 * <code>propto</code> parameter is <code>false</code> or if any
 * of the other template arguments are not constants as defined by
 * <code>stan::is_constant_all<T></code>.
 *
 * Example use: <code>include_summand<false, double, var, double, double></code>
 *
 * @tparam propto <code>true</code> if calculating up to a
 * proportionality constant.
 * @tparam T (optional). A type
 * @tparam T_pack (optional). A parameter pack of types. This is used to
 * extend the applicability of the function to an arbitrary number of types.
 */
template <bool propto, typename T = double, typename... T_pack>
struct include_summand
    : bool_constant<(!stan::is_constant_all<scalar_type_t<T>>::value
                     || include_summand<propto, T_pack...>::value)> {};

/** \ingroup type_trait
 * <code>true</code> if a term with the specified propto
 * value and subterm types should be included in a proportionality
 * calculation.
 */
template <bool propto, typename T>
struct include_summand<propto, T>
    : bool_constant<(
          !propto
          || !stan::is_constant_all<typename scalar_type<T>::type>::value)> {};

}  // namespace math

}  // namespace stan

#endif
