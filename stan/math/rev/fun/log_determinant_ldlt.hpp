#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_LDLT_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the log det of the matrix whose LDLT factorization is given
 *
 * @param A an LDLT_factor2
 * @return ln(det(A))
 */
template <typename T, bool alloc_in_arena,
	  require_rev_matrix_t<T>* = nullptr>
var log_determinant_ldlt(const LDLT_factor2<T, alloc_in_arena> &A) {
  if (A.matrix().size() == 0) {
    return 0;
  }
  
  var log_det = sum(log(A.ldlt().vectorD().array()));

  arena_t<Eigen::MatrixXd> arena_A_inv(A.matrix().rows(),
				       A.matrix().cols());

  arena_A_inv.setIdentity();
  A.ldlt().solveInPlace(arena_A_inv);

  reverse_pass_callback([A, log_det, arena_A_inv]() mutable {
    A.matrix().adj() += log_det.adj() * arena_A_inv;
  });

  return log_det;
}
  
}  // namespace math
}  // namespace stan
#endif
