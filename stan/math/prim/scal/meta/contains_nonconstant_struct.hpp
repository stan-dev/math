#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_NONCONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_NONCONSTANT_STRUCT_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/or.hpp>
#include <type_traits>
namespace stan {
/**
 * Metaprogram to determine if any of the
 * provided types have a base scalar
 * type that can not be assigned to type double.
 */
template <typename... T>
using contains_nonconstant_struct = math::or_not_<is_constant_struct<T>...>;

}  // namespace stan
#endif
