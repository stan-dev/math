#ifndef STAN_MATH_PRIM_FUN_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_IDENTITY_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * <p>This method is effectively a no-op and is mainly useful as a
 * placeholder in auto-generated code.
 *
 * @tparam T Any type.
 * @param[in] x object
 * @return transformed input
 */
template <typename T>
inline auto identity_constrain(T&& x) {
  return std::forward<T>(x);
}

/**
 * Returns the result of applying the identity constraint
 * transform to the input and increments the log probability
 * reference with the log absolute Jacobian determinant.
 *
 * <p>This method is effectively a no-op and mainly useful as a
 * placeholder in auto-generated code.
 *
 * @tparam T Any type
 * @tparam S type of log probability
 * @param[in] x object
 * @param[in] lp log density reference
 * @return transformed input
 */
template <typename T, typename S>
inline auto identity_constrain(T&& x, S& lp) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan

#endif
