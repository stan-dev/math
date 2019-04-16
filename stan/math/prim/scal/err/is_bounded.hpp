#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_BOUNDED_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_BOUNDED_HPP

#include <stan/math/prim/scal/meta/max_size.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/get.hpp>

namespace stan {
namespace math {

/**
 * Check if the value is between the low and high values, inclusively.
 * @tparam T_y Type of value
 * @tparam T_low Type of low value
 * @tparam T_high Type of high value
 * @param y Value to check
 * @param low Low bound
 * @param high High bound
 * @return <code>true</code> if the value provided is within the bound provided
 *   and none of the arguments are NaN
 */
template <typename T_y, typename T_low, typename T_high>
inline bool is_bounded(const T_y& y, const T_low& low, const T_high& high) {
  scalar_seq_view<T_low> low_vec(low);
  scalar_seq_view<T_high> high_vec(high);
  for (size_t n = 0; n < stan::max_size(low, high); n++) {
    if (!(low_vec[n] <= y && y <= high_vec[n]))
      return false;
  }

  for (size_t n = 0; n < stan::length(y); n++) {
    if (!(low_vec[n] <= stan::get(y, n) && stan::get(y, n) <= high_vec[n]))
      return false;
  }

  return true;
}

}  // namespace math
}  // namespace stan
#endif
