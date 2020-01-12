#ifndef STAN_MATH_PRIM_ERR_IS_SYMMETRIC_HPP
#define STAN_MATH_PRIM_ERR_IS_SYMMETRIC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/is_square.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {
/**
 * Return <code>true</code> if the matrix is square, and no element
 * not on the main diagonal is <code>NaN</code>.
 * @tparam T_y Type of scalar, requires class method <code>.rows()</code>
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square, and no
 *    element not on the main diagonal is <code>NaN</code>
 */
template <typename T_y>
inline bool is_symmetric(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  if (!is_square(y)) {
    return false;
  }

  using size_type = typename index_type<
      Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>>::type;

  size_type k = y.rows();
  if (k == 1) {
    return true;
  }
  for (size_type m = 0; m < k; ++m) {
    for (size_type n = m + 1; n < k; ++n) {
      if (!(fabs(value_of(y(m, n)) - value_of(y(n, m)))
            <= CONSTRAINT_TOLERANCE)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
