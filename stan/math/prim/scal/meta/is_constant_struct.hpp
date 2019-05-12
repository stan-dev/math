#ifndef STAN_MATH_PRIM_SCAL_META_IS_CONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_CONSTANT_STRUCT_HPP

#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>

namespace stan {

/**
 * Metaprogram to determine if a type has a base scalar
 * type that can be assigned to type double.
 * @tparam T Types to test
 */
template <typename... T>
using is_constant_struct = math::conjunction<is_constant<T>...>;

}  // namespace stan
#endif
