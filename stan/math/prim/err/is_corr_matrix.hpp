#ifndef STAN_MATH_PRIM_ERR_IS_CORR_MATRIX_HPP
#define STAN_MATH_PRIM_ERR_IS_CORR_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/is_pos_definite.hpp>
#include <stan/math/prim/err/is_positive.hpp>
#include <stan/math/prim/err/is_size_match.hpp>
#include <stan/math/prim/err/is_symmetric.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the matrix is square and not 0x0,
 * if the matrix is symmetric, diagonals are near 1, positive definite,
 * and no elements are <code>NaN</code>
 * A valid correlation matrix is symmetric, has a unit diagonal
 * (all 1 values), and has all values between -1 and 1 (inclusive).
 * @tparam T_y Type of scalar, requires class method <code>.rows()</code>
 *   and <code>.cols()</code>
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square and not 0x0,
 *   if the matrix is symmetric, diagonals are near 1, positive definite,
 *   and no elements are <code>NaN</code>
 */
template <typename T_y>
inline bool is_corr_matrix(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  using size_type
      = index_type_t<Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>>;
  using std::fabs;

  if (!is_size_match(y.rows(), y.cols())) {
    return false;
  }
  if (!is_positive(y.rows())) {
    return false;
  }
  if (!is_pos_definite(y)) {
    return false;
  }
  if (is_symmetric(y)) {
    for (size_type k = 0; k < y.rows(); ++k) {
      if (!(fabs(y(k, k) - 1.0) <= CONSTRAINT_TOLERANCE)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
