#ifndef STAN_MATH_PRIM_SCAL_FUN_PROMOTE_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_PROMOTE_SCALAR_TYPE_HPP

#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * Template metaprogram to calculate a type for converting a
 * convertible type.  This is the base case.
 *
 * @tparam T result scalar type.
 * @tparam S input type
 */
template <typename T, typename S>
struct promote_scalar_type {
  /**
   * The promoted type.
   */
  using type = T;
};

template <typename T, typename S>
using promote_scalar_t = typename promote_scalar_type<T, S>::type;

}  // namespace math
}  // namespace stan
#endif
