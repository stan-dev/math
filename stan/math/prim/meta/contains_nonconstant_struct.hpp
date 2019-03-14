#ifndef STAN_MATH_PRIM_META_CONTAINS_NONCONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_META_CONTAINS_NONCONSTANT_STRUCT_HPP

#include <stan/math/prim/meta/is_nonconstant_struct.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Metaprogram to determine if any of the
 * provided types have a base scalar
 * type that cannot be assigned to type double.
 */
template <typename... T>
using contains_nonconstant_struct
    = math::disjunction<is_nonconstant_struct<T>...>;

}  // namespace stan
#endif
