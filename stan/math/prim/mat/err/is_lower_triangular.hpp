#ifndef STAN_MATH_PRIM_MAT_ERR_IS_LOWER_TRIANGULAR_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_LOWER_TRIANGULAR_HPP

#include <Eigen/Dense>
#include <math.h>

namespace stan {
namespace math {

double notLowerNan(double x) { return std::isnan(x) ? 1.0 : x; }

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
  Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic> m
      = y.unaryExpr(std::ptr_fun(notLowerNan));
  return m.transpose().isUpperTriangular();
}

}  // namespace math
}  // namespace stan
#endif
