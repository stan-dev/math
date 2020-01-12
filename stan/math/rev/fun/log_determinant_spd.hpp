#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
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
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param m a symmetric, positive-definite matrix
 * @return The log determinant of the specified matrix
 */
template <int R, int C>
inline var log_determinant_spd(const Eigen::Matrix<var, R, C>& m) {
  check_symmetric("log_determinant_spd", "m", m);
  if (m.size() == 0)
    return 0;

  matrix_d m_d = m.val();

  Eigen::LDLT<matrix_d> ldlt(m_d);
  if (ldlt.info() != Eigen::Success) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "failed LDLT factorization");
  }

  // compute the inverse of A (needed for the derivative)
  m_d.setIdentity(m.rows(), m.cols());
  ldlt.solveInPlace(m_d);

  if (ldlt.isNegative() || (ldlt.vectorD().array() <= 1e-16).any()) {
    double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "matrix is negative definite");
  }

  double val = sum(log(ldlt.vectorD()));

  check_finite("log_determinant_spd",
               "log determininant of the matrix argument", val);

  vari** operands
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(m.size());
  Eigen::Map<matrix_vi>(operands, m.rows(), m.cols()) = m.vi();

  double* gradients
      = ChainableStack::instance_->memalloc_.alloc_array<double>(m.size());
  Eigen::Map<matrix_d>(gradients, m.rows(), m.cols()) = m_d;

  return var(
      new precomputed_gradients_vari(val, m.size(), operands, gradients));
}

}  // namespace math
}  // namespace stan
#endif
