#ifndef STAN_MATH_PRIM_META_CONTAINS_VECTOR_HPP
#define STAN_MATH_PRIM_META_CONTAINS_VECTOR_HPP

#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/disjunction.hpp>

namespace stan {
/** \ingroup type_trait
 * Metaprogram to determine if any of the
 * provided types is a std or eigen vector.
 * @tparam T Types to test
 */
template <typename... T>
using contains_vector = math::disjunction<is_vector<T>...>;

}  // namespace stan
#endif
