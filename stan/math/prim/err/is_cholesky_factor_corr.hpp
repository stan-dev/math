#ifndef STAN_MATH_PRIM_ERR_IS_CHOLESKY_FACTOR_CORR_HPP
#define STAN_MATH_PRIM_ERR_IS_CHOLESKY_FACTOR_CORR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/err/is_cholesky_factor.hpp>
#include <stan/math/prim/err/is_unit_vector.hpp>

namespace stan {
namespace math {
/**
 * Return <code>true</code> if y is a valid Cholesky factor, if
 * the number of rows is not less than the number of columns, if there
 * are no zero columns, and no element in matrix is <code>NaN</code>.
 * A Cholesky factor is a lower triangular matrix whose diagonal
 * elements are all positive. This definition does not require a
 * square matrix.
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param y Matrix to test
 * @return <code>true</code> if y is a valid Cholesky factor, if
 *    the number of rows is not less than the number of columns,
 *    if there are no 0 columns, and no element in matrix is <code>NaN</code>
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline bool is_cholesky_factor_corr(const EigMat& y) {
  const auto& y_ref = to_ref(y);
  if (!is_cholesky_factor(y_ref)) {
    return false;
  }
  for (int i = 0; i < y_ref.rows(); ++i) {
    if (!is_unit_vector(y_ref.row(i))) {
      return false;
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
