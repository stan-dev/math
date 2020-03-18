#ifndef STAN_MATH_PRIM_META_IS_CONTAINER_HPP
#define STAN_MATH_PRIM_META_IS_CONTAINER_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Deduces whether type is eigen matrix or standard vector.
 * @tparam Container type to check
 */
template <typename Container>
using is_container = bool_constant<
    math::disjunction<is_eigen<Container>, is_std_vector<Container>>::value>;

/** \addtogroup require_container_types
*  @{
*/
/**
 * Requires a type is an eigen or std vector.
 */
STAN_ADD_REQUIRE_UNARY(is_container, is_container);
STAN_ADD_REQUIRE_UNARY_SCALAR(is_container, is_container);
STAN_ADD_REQUIRE_UNARY_VALUE(is_container, is_container);
/** @}*/

}  // namespace stan

#endif
