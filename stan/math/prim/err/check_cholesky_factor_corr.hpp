#ifndef STAN_MATH_PRIM_ERR_CHECK_CHOLESKY_FACTOR_CORR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CHOLESKY_FACTOR_CORR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/err/check_lower_triangular.hpp>
#include <stan/math/prim/err/check_square.hpp>
#include <stan/math/prim/err/check_unit_vector.hpp>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is a valid Cholesky factor of a
 * correlation matrix.
 * A Cholesky factor is a lower triangular matrix whose diagonal
 * elements are all positive.  Note that Cholesky factors need not
 * be square, but require at least as many rows M as columns N
 * (i.e., M &gt;= N).
 * Tolerance is specified by <code>math::CONSTRAINT_TOLERANCE</code>.
 * @tparam EigMat Type inheriting from `MatrixBase` with dynamic rows and
 * columns.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw <code>std::domain_error</code> if y is not a valid Cholesky
 *   factor, if number of rows is less than the number of columns,
 *   if there are 0 columns, or if any element in matrix is NaN
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
void check_cholesky_factor_corr(const char* function, const char* name,
                                const EigMat& y) {
  const auto& y_ref = to_ref(y);
  check_square(function, name, y_ref);
  check_lower_triangular(function, name, y_ref);
  check_positive(function, name, y_ref.diagonal());
  for (Eigen::Index i = 0; i < y_ref.rows(); ++i) {
    check_unit_vector(function, name, y_ref.row(i));
  }
}

}  // namespace math
}  // namespace stan
#endif
