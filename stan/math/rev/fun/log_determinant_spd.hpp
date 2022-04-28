#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
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
 * @tparam T Type is an Eigen Matrix
 * @param m a symmetric, positive-definite matrix
 * @return The log determinant of the specified matrix
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var log_determinant_spd(const T& m) {
  if (m.size() == 0) {
    return var(0.0);
  }
  check_symmetric("log_determinant", "m", m);

  matrix_d m_d = m.val();
  arena_t<T> arena_m = m;
  auto m_ldlt = arena_m.val().ldlt();
  if (m_ldlt.info() != Eigen::Success) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "failed LDLT factorization");
  }

  // compute the inverse of A (needed for the derivative)
  m_d.setIdentity(m.rows(), m.cols());

  auto arena_m_inv_transpose = to_arena(m_ldlt.solve(m_d).transpose());

  if (m_ldlt.isNegative() || (m_ldlt.vectorD().array() <= 1e-16).any()) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "matrix is negative definite");
  }

  var log_det = sum(log(m_ldlt.vectorD()));

  reverse_pass_callback([arena_m, log_det, arena_m_inv_transpose]() mutable {
    arena_m.adj() += log_det.adj() * arena_m_inv_transpose;
  });
  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
