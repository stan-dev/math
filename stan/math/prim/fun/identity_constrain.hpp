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
 * @tparam T Type of scalar.
 * @param[in] x free scalar
 * @return transformed input
 */
template <typename T>
inline T identity_constrain(const T& x) {
  return x;
}

/**
 * Returns the result of applying the identity constraint
 * transform to the input and increments the log probability
 * reference with the log absolute Jacobian determinant.
 *
 * <p>This method is effectively a no-op and mainly useful as a
 * placeholder in auto-generated code.
 *
 * @tparam T type of scalar
 * @tparam S type of log probability
 * @param[in] x scalar
 * @param[in] lp log density reference
 * @return transformed input
 */
template <typename T, typename S>
inline T identity_constrain(const T& x, S& lp) {
  return x;
}

}  // namespace math
}  // namespace stan

#endif
