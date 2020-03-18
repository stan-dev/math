#ifndef STAN_MATH_PRIM_META_IS_DOUBLE_OR_INT_HPP
#define STAN_MATH_PRIM_META_IS_DOUBLE_OR_INT_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {
  /**
   * Checks if decayed type is a double or integer
   * @tparam The type to check
   */
  template <typename T>
  struct is_double_or_int
      : bool_constant<
            math::disjunction<std::is_same<double, std::decay_t<T>>,
                              std::is_same<int, std::decay_t<T>>>::value> {};

/** \addtogroup require_stan_scalar
*  @{
*/
/**
 * Require that a type is `double` or `int`.
 */
STAN_ADD_REQUIRE_UNARY(double_or_int, is_double_or_int);
STAN_ADD_REQUIRE_UNARY_SCALAR(double_or_int, is_double_or_int);
STAN_ADD_REQUIRE_UNARY_VALUE(double_or_int, is_double_or_int);
/** @}*/

}  // namespace stan
#endif
