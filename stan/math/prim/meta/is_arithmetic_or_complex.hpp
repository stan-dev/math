#ifndef STAN_MATH_PRIM_META_IS_ARITHMETIC_OR_COMPLEX_HPP
#define STAN_MATH_PRIM_META_IS_ARITHMETIC_OR_COMPLEX_HPP

#include <stan/math/prim/meta/is_complex.hpp>
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
struct is_arithmetic_or_complex
    : math::disjunction<std::is_arithmetic<std::decay_t<T>>,
                                      is_complex<std::decay_t<T>>> {};

STAN_ADD_REQUIRE_UNARY(arithmetic_or_complex, is_arithmetic_or_complex, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(arithmetic_or_complex, is_arithmetic_or_complex, require_stan_scalar_real);

}  // namespace stan

#endif
