#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
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

  arena_t<T> arena_m = m;

  auto arena_m_inv = to_arena(arena_m.val());
  check_symmetric("log_determinant_spd", "m", arena_m_inv);

  const auto& ldlt = arena_m_inv.ldlt();

  if (ldlt.info() != Eigen::Success) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "failed LDLT factorization");
  }

  if (ldlt.isNegative() || (ldlt.vectorD().array() <= 1e-16).any()) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "matrix is negative definite");
  }

  // compute the inverse of A (needed for the derivative)
  arena_m_inv.setIdentity(m.rows(), m.cols());
  ldlt.solveInPlace(arena_m_inv);

  var log_det = sum(log(ldlt.vectorD()));

  check_finite("log_determinant_spd",
               "log determininant of the matrix argument", log_det);

  reverse_pass_callback([arena_m, log_det, arena_m_inv]() mutable {
    arena_m.adj() += log_det.adj() * arena_m_inv;
  });

  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
