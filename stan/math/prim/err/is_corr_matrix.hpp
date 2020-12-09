#ifndef STAN_MATH_PRIM_ERR_IS_CORR_MATRIX_HPP
#define STAN_MATH_PRIM_ERR_IS_CORR_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/is_pos_definite.hpp>
#include <stan/math/prim/err/is_positive.hpp>
#include <stan/math/prim/err/is_size_match.hpp>
#include <stan/math/prim/err/is_symmetric.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the matrix is square and not 0x0,
 * if the matrix is symmetric, diagonals are near 1, positive definite,
 * and no elements are <code>NaN</code>
 * A valid correlation matrix is symmetric, has a unit diagonal
 * (all 1 values), and has all values between -1 and 1 (inclusive).
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square and not 0x0,
 *   if the matrix is symmetric, diagonals are near 1, positive definite,
 *   and no elements are <code>NaN</code>
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline bool is_corr_matrix(const EigMat& y) {
  using std::fabs;
  const auto& y_ref = to_ref(y);
  if (!is_size_match(y_ref.rows(), y_ref.cols())) {
    return false;
  }
  if (!is_positive(y_ref.rows())) {
    return false;
  }
  if (!is_pos_definite(y_ref)) {
    return false;
  }
  if (is_symmetric(y_ref)) {
    for (Eigen::Index k = 0; k < y_ref.rows(); ++k) {
      if (!(fabs(y_ref(k, k) - 1.0) <= CONSTRAINT_TOLERANCE)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
