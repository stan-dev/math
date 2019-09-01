#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

namespace stan {
/**
 * Metaprogram to determine if any of the
 * provided types is a std or eigen vector.
 * @tparam T Types to test
 */
template <typename... T>
using contains_vector = math::disjunction<is_vector<T>...>;

template<class... T>
constexpr bool is_contains_vector_v = contains_vector<T...>::value;

}  // namespace stan
#endif
