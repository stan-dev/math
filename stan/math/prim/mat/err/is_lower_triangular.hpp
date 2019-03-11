#ifndef STAN_MATH_PRIM_MAT_ERR_IS_LOWER_TRIANGULAR_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_LOWER_TRIANGULAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> is matrix is lower triangular.
 * A matrix x is not lower triangular if there is a non-zero entry
 * x[m, n] with m &lt; n. This function only inspect the upper and
 * triangular portion of the matrix, not including the diagonal.
 * @tparam T Type of scalar of the matrix, requires class method
 *   <code>.rows()</code> and <code>.cols()</code>
 * @param y Matrix to test
 * @return <code>true</code> is matrix is lower triangular
 */
template <typename T_y>
inline bool is_lower_triangular(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  for (int n = 1; n < y.cols(); ++n) {
    for (int m = 0; m < n && m < y.rows(); ++m) {
      if (y(m, n) != 0)
        return false;
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
