#ifndef STAN_MATH_PRIM_MAT_ERR_IS_POSITIVE_ORDERED_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_POSITIVE_ORDERED_HPP

#include <stan/math/prim/mat/err/is_ordered.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified vector contains only non-negative
 * values and is sorted into strictly increasing order.
 * @tparam T_y Type of the variable, requires class method <code>.size()</code>
 * @param y Vector to test
 * @return <code>true</code> if the vector contains only non-negative
 *   values, if the values are ordered, if there are no duplicated
 *   values, and if no element is <code>NaN</code>
 */
template <typename T_y>
inline bool is_positive_ordered(
    const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y) {
  if (y.size() == 0)
    return true;
  if (y[0] < 0)
    return false;
  return is_ordered(y);
}

}  // namespace math
}  // namespace stan
#endif
