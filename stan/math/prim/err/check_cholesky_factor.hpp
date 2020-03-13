#ifndef STAN_MATH_PRIM_ERR_CHECK_CHOLESKY_FACTOR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CHOLESKY_FACTOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/err/check_less_or_equal.hpp>
#include <stan/math/prim/err/check_lower_triangular.hpp>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is a valid Cholesky factor.
 * A Cholesky factor is a lower triangular matrix whose diagonal
 * elements are all positive.  Note that Cholesky factors need not
 * be square, but require at least as many rows M as columns N
 * (i.e., M &gt;= N).
 * @tparam EigMat Type of the Cholesky factor (must be derived from \c
 * Eigen::MatrixBase)
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw <code>std::domain_error</code> if y is not a valid Cholesky
 *   factor, if number of rows is less than the number of columns,
 *   if there are 0 columns, or if any element in matrix is NaN
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline void check_cholesky_factor(const char* function, const char* name,
                                  const EigMat& y) {
  check_less_or_equal(function, "columns and rows of Cholesky factor", y.cols(),
                      y.rows());
  check_positive(function, "columns of Cholesky factor", y.cols());
  const Eigen::Ref<const plain_type_t<EigMat>>& y_ref = y;
  check_lower_triangular(function, name, y_ref);
  check_positive(function, name, y_ref.diagonal());
}

}  // namespace math
}  // namespace stan
#endif
