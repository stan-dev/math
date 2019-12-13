#ifndef STAN_MATH_PRIM_META_CONTAINS_FVAR_HPP
#define STAN_MATH_PRIM_META_CONTAINS_FVAR_HPP

#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/disjunction.hpp>

namespace stan {

/** \ingroup type_trait
 * Extends std::true_type when instantiated with at least 1
 * template parameter that is a fvar. Extends std::false_type
 * otherwise.
 * @tparam T Types to test
 */
template <typename... T>
using contains_fvar = math::disjunction<is_fvar<scalar_type_t<T>>...>;

}  // namespace stan
#endif
