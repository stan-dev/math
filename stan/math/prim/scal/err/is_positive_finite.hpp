#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_POSITIVE_FINITE_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_POSITIVE_FINITE_HPP

#include <stan/math/prim/scal/err/is_positive.hpp>
#include <stan/math/prim/scal/err/is_scal_finite.hpp>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is positive and finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of <code>y</code>
 * @param y Variable to check
 * @return <code>true</code> if any element of <code>y</code> is positive and
 *   if no element of <code>y</code> is NaN
 */
template <typename T_y>
inline bool is_positive_finite(const T_y& y) {
  return is_positive(y) && is_scal_finite(y);
}

}  // namespace math
}  // namespace stan
#endif
