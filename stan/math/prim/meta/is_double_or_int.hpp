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
 * @ingroup type_trait
 */
template <typename T>
struct is_double_or_int
    : bool_constant<
          math::disjunction<std::is_same<double, std::decay_t<T>>,
                            std::is_same<int, std::decay_t<T>>>::value> {};

STAN_ADD_REQUIRE_UNARY(double_or_int, is_double_or_int,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(double_or_int, is_double_or_int,
                             require_stan_scalar_real);

}  // namespace stan
#endif
