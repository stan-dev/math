#ifndef STAN_MATH_PRIM_MAT_ERR_IS_CHOLESKY_FACTOR_CORR_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_CHOLESKY_FACTOR_CORR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/is_positive.hpp>
#include <stan/math/prim/mat/err/is_lower_triangular.hpp>
#include <stan/math/prim/mat/err/is_square.hpp>
#include <stan/math/prim/mat/err/constraint_tolerance.hpp>
#include <stan/math/prim/mat/err/is_unit_vector.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if y is a valid Cholesky factor, if
 * the number of rows is not less than the number of columns, if there
 * are no zero columns, and no element in matrix is <code>NaN</code>.
 * A Cholesky factor is a lower triangular matrix whose diagonal
 * elements are all positive. This definition does not require a
 * square matrix just that M &gt;= N, for M rows and N columns.
 * Tolerance is specified by <code>math::CONSTRAINT_TOLERANCE</code>
 * as 1E-8.
 * @tparam T_y Type of elements of Cholesky factor, requires class method
 *   <code>.rows()</code>
 * @param y Matrix to test
 * @return <code>true</code> if y is a valid Cholesky factor, if
 *    the number of rows is not less than the number of columns,
 *    if there are no 0 columns, and no element in matrix is <code>NaN</code>
 */
template <typename T_y>
inline bool is_cholesky_factor_corr(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  if (is_square(y)) {
    if (is_lower_triangular(y)) {
      for (int i = 0; i < y.rows(); ++i) {
        if (!is_positive(y(i, 1)))
          return false;
      }
      for (int i = 0; i < y.rows(); ++i) {
        Eigen::Matrix<T_y, Eigen::Dynamic, 1> y_i = y.row(i).transpose();
        if (!is_unit_vector(y_i))
          return false;
      }
      return true;
    }
  }
  return false;
}

}  // namespace math
}  // namespace stan
#endif
