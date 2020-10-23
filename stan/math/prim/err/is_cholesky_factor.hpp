#ifndef STAN_MATH_PRIM_ERR_IS_CHOLESKY_FACTOR_HPP
#define STAN_MATH_PRIM_ERR_IS_CHOLESKY_FACTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_positive.hpp>
#include <stan/math/prim/err/is_lower_triangular.hpp>
#include <stan/math/prim/err/is_less_or_equal.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if y is a valid Cholesky factor, if
 * number of rows is not less than the number of columns, if there
 * are no 0 columns, and no element in matrix is <code>NaN</code>.
 * A Cholesky factor is a lower triangular matrix whose diagonal
 * elements are all positive.  Note that Cholesky factors need not
 * be square, but require at least as many rows M as columns N
 * (i.e., M &gt;= N).
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param y Matrix to test
 * @return <code>true</code> if y is a valid Cholesky factor, if
 *   number of rows is not less than the number of columns,
 *   if there are no 0 columns, and no element in matrix is <code>NaN</code>
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline bool is_cholesky_factor(const EigMat& y) {
  const auto& y_ref = to_ref(y);
  return is_less_or_equal(y_ref.cols(), y_ref.rows())
         && is_positive(y_ref.cols()) && is_lower_triangular(y_ref)
         && is_positive(y_ref.diagonal());
}

}  // namespace math
}  // namespace stan
#endif
