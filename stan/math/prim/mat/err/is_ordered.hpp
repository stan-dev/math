#ifndef STAN_MATH_PRIM_MAT_ERR_IS_ORDERED_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_ORDERED_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>
#include <stan/math/prim/scal/meta/error_index.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified vector is sorted into
 * strictly increasing order.
 * @tparam T_y Type of scalar, require class method <code>.size()</code>
 * @param y Vector to test
 * @return <code>true</code> if the vector elements are ordered, if 
 *   there are no duplicated values, and if no element is <code>NaN</code>
 */
template <typename T_y>
inline bool is_ordered(const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  typedef typename index_type<Matrix<T_y, Dynamic, 1> >::type size_t;

  for (size_t n = 1; n < y.size(); ++n) {
    if (!(y[n] > y[n - 1])) {
      return false;
    }
  } return true;
}

}  // namespace math
}  // namespace stan
#endif
