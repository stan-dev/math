#ifndef STAN_MATH_PRIM_META_IS_STAN_SCALAR_HPP
#define STAN_MATH_PRIM_META_IS_STAN_SCALAR_HPP

#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is a var, fvar, or arithmetic
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_stan_scalar
    : bool_constant<
          math::disjunction<is_var<std::decay_t<T>>, is_fvar<std::decay_t<T>>,
                            std::is_arithmetic<std::decay_t<T>>>::value> {};

STAN_ADD_REQUIRE_UNARY(stan_scalar, is_stan_scalar, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(stan_scalar, is_stan_scalar,
                             require_stan_scalar_real);

}  // namespace stan

#endif
