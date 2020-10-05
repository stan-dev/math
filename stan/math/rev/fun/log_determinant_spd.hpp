#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Returns the log det of a symmetric, positive-definite matrix
 *
 * @tparam EigMat Type of the matrix
 * @param m a symmetric, positive-definite matrix
 * @return The log determinant of the specified matrix
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline var log_determinant_spd(const T& m) {
  if (m.size() == 0) {
    return 0;
  }

  arena_matrix<Eigen::MatrixXd> arena_m_d = m.val();
  check_symmetric("log_determinant_spd", "m", arena_m_d);

  Eigen::LDLT<Eigen::MatrixXd> ldlt(arena_m_d);
  if (ldlt.info() != Eigen::Success) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "failed LDLT factorization");
  }

  // compute the inverse of A (needed for the derivative)
  arena_m_d.setIdentity(m.rows(), m.cols());
  ldlt.solveInPlace(arena_m_d);

  if (ldlt.isNegative() || (ldlt.vectorD().array() <= 1e-16).any()) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "matrix is negative definite");
  }

  var log_det = sum(log(ldlt.vectorD()));

  check_finite("log_determinant_spd",
               "log determininant of the matrix argument", log_det);

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arena_m = m;

  reverse_pass_callback([arena_m, log_det, arena_m_d]() mutable {
    arena_m.adj() += log_det.adj() * arena_m_d;
  });

  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
