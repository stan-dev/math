#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_LESS_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_LESS_HPP

#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/err/domain_error_vec.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/get.hpp>
#include <functional>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is strictly less than <code>high</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>high</code>.
 * @tparam T_y Type of <code>y</code>
 * @tparam T_high Type of upper bound
 * @param y Variable to check
 * @param high Upper bound
 * @return <code>true</code> if <code>y</code> is less than low and no element
 *   of <code>y</code> or <code>high</code> is NaN
 */
template <typename T_y, typename T_high>
inline bool is_less(const T_y& y, const T_high& high) {
  scalar_seq_view<T_high> high_vec(high);
  for (size_t n = 0; n < stan::length(high); n++) {
    if (!(y < high_vec[n]))
      return false;
  }

  for (size_t n = 0; n < stan::length(y); n++) {
    if (!(stan::get(y, n) < high_vec[n]))
      return false;
  }

  return true;
}

}  // namespace math
}  // namespace stan
#endif
