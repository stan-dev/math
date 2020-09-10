#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_LDLT_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the log det of the matrix whose LDLT factorization is given
 * See The Matrix Cookbook's chapter on Derivatives of a Determinant
 * In this case, it is just the inverse of the underlying matrix
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param A an LDLT_factor
 * @return ln(det(A))
 */
template <int R, int C>
var log_determinant_ldlt(LDLT_factor<var, R, C> &A) {
  if (A.rows() == 0) {
    return 0;
  }

  var log_det = A.alloc_->log_abs_det();

  arena_matrix<Eigen::Matrix<double, R, C>> arena_A_inv(A.rows(), A.cols());

  arena_A_inv.setIdentity();
  A.alloc_->ldlt_.solveInPlace(arena_A_inv);

  reverse_pass_callback([A, log_det, arena_A_inv]() mutable {
    A.alloc_->arena_A_.adj() += log_det.adj() * arena_A_inv;
  });

  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
