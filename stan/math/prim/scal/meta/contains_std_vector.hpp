#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_STD_VECTOR_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <type_traits>

namespace stan {
/**
 * Checks if any types are std vectors.
 */
template <typename... Ts>
using contains_std_vector = math::disjunction<is_std_vector<Ts>...>;
}  // namespace stan

#endif
