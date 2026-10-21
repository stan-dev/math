#ifndef STAN_MATH_PRIM_META_CONTAINS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_META_CONTAINS_STD_VECTOR_HPP

#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_vector.hpp>

namespace stan {
/** \ingroup type_trait
 * Checks if any types are std vectors.
 */
template <typename... Ts>
using contains_std_vector = math::disjunction<is_std_vector<Ts>...>;
}  // namespace stan

#endif
