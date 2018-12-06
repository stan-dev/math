#ifndef STAN_MATH_PRIM_MAT_ERR_CHECK_POS_SEMIDEFINITE_HPP
#define STAN_MATH_PRIM_MAT_ERR_CHECK_POS_SEMIDEFINITE_HPP

#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/mat/err/constraint_tolerance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <sstream>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is positive definite
 *
 * @tparam T_y scalar type of the matrix
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 *
 * @throw <code>std::invalid_argument</code> if the matrix is not square
 *   or if the matrix has 0 size.
 * @throw <code>std::domain_error</code> if the matrix is not symmetric,
 *   or if it is not positive semi-definite,
 *   or if any element of the matrix is <code>NaN</code>.
 */
template <typename T_y>
inline void check_pos_semidefinite(
    const char* function, const char* name,
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  check_symmetric(function, name, y);
  check_positive_size(function, name, "rows", y.rows());

  if (y.rows() == 1 && !(y(0, 0) >= 0.0))
    domain_error(function, name, "is not positive semi-definite.", "");

  using Eigen::Dynamic;
  using Eigen::LDLT;
  using Eigen::Matrix;
  LDLT<Matrix<double, Dynamic, Dynamic> > cholesky = value_of_rec(y).ldlt();
  if (cholesky.info() != Eigen::Success
      || (cholesky.vectorD().array() < 0.0).any())
    domain_error(function, name, "is not positive semi-definite.", "");
  check_not_nan(function, name, y);
}

/**
 * Check if the specified matrix is positive semidefinite
 *
 * @tparam Derived Derived type of the Eigen::LDLT transform.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param cholesky Eigen::LDLT to test, whose progenitor matrix must
 * not have any NaN elements.
 * @throw <code>std::invalid_argument</code> if the matrix has 0 size.
 * @throw <code>std::domain_error</code> if the matrix is not positive
 *   semi-definite.
 */
template <typename Derived>
inline void check_pos_semidefinite(
    const char* function, const char* name,
    const Eigen::LDLT<Derived>& cholesky) {
  // From the Eigen::LDLT we cannot check for NaNs in the original
  // matrix, nor can we test for symmetry. Eigen::LDLT assumes the
  // matrix is symmetric and uses only half of it.
  // We also cannot check if that matrix was size 0. Eigen::LDLT on a
  // size zero matrix can give a segfault, although this error
  // behavior does not seem to be documented. So it is on the caller
  // to handle this case.
  if (cholesky.info() != Eigen::Success
      || (cholesky.vectorD().array() < 0.0).any())
    domain_error(function, name, "is not positive semi-definite.", "");
}

}  // namespace math
}  // namespace stan
#endif
