#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_GREATER_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_GREATER_HPP

#include <stan/math/prim/meta.hpp>
#include <functional>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is strictly greater than <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y Type of <code>y</code>
 * @tparam T_low Type of lower bound
 * @param y Variable to check
 * @param low Lower bound
 * @return <code>true</code> if <code>y</code> is greater than low and if no
 *   element of <code>y</code> or <code>low</code> is NaN.
 */
template <typename T_y, typename T_low>
inline bool is_greater(const T_y& y, const T_low& low) {
  scalar_seq_view<T_low> low_vec(low);
  for (size_t n = 0; n < stan::length(low); n++) {
    if (!(y > low_vec[n]))
      return false;
  }

  for (size_t n = 0; n < stan::length(y); n++) {
    if (!(stan::get(y, n) > low_vec[n]))
      return false;
  }

  return true;
}

}  // namespace math
}  // namespace stan
#endif
