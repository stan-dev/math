#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Returns the scalar type inside of a possibly nested container.
 * - Scalar types this return a scalar such as `double`, `int` or `var`
 * - Eigen types this will return their `Scalar` member
 * - `std::vector<T>` will call `scalar_type<T>` which will
 *  continue to be called recursively until a scalar type is reached.
 * The recursive part of `scalar_type` is what seperates it from `value_type`
 * in that `value_type` will only go one level down in a container of
 * containers while `scalar_type` will recursively call `scalar_type` on
 * inner types until a scalar value is reached.
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T type of non-container
 */
template <typename T, typename = void>
struct scalar_type {
  using type = std::decay_t<T>;
};

/** \ingroup type_trait
 * Helper function. See the docs for `scalar_type`
 */
template <typename T>
using scalar_type_t = typename scalar_type<T>::type;

}  // namespace stan
#endif
