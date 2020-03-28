#ifndef STAN_MATH_PRIM_META_IS_AUTODIFF_HPP
#define STAN_MATH_PRIM_META_IS_AUTODIFF_HPP

#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <complex>
#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_autodiff
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

STAN_ADD_REQUIRE_UNARY(autodiff, is_autodiff, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(autodiff, is_autodiff, require_stan_scalar_real);

}  // namespace stan

#endif
