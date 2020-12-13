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
 * @tparam EigMat Type of the matrix
 * @param m a symmetric, positive-definite matrix
 * @return The log determinant of the specified matrix
 */
template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline var log_determinant_spd(const EigMat& m) {
  if (m.size() == 0) {
    return 0;
  }

  matrix_d m_d = m.val();
  check_symmetric("log_determinant_spd", "m", m_d);

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

/**
 * Returns the log det of a symmetric, positive-definite matrix
 *
 * @tparam EigMat Type of the matrix
 * @param m a symmetric, positive-definite matrix
 * @return The log determinant of the specified matrix
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var log_determinant_spd(const T& m) {
  check_square("log_determinant_spd", "m", m);

  if (m.size() == 0) {
    return var(0.0);
  }

  check_symmetric("log_determinant_spd", "m", m.val());

  auto ldlt = m.val().ldlt();
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

  arena_t<Eigen::MatrixXd> arena_m_inv_transpose
      = Eigen::MatrixXd::Identity(m.rows(), m.cols());
  ldlt.solveInPlace(arena_m_inv_transpose);

  return make_callback_var(sum(log(ldlt.vectorD())),
                           [m, arena_m_inv_transpose](const auto& res) mutable {
                             m.adj() += res.adj() * arena_m_inv_transpose;
                           });
}

}  // namespace math
}  // namespace stan
#endif
