#ifndef STAN_MATH_REV_MAT_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_REV_MAT_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

template <int R, int C>
inline var log_determinant_spd(const Eigen::Matrix<var, R, C>& m) {
  check_square("log_determinant_spd", "m", m);

  matrix_d m_d = m.val();

  Eigen::LDLT<matrix_d> ldlt(m_d);
  if (ldlt.info() != Eigen::Success) {
    double y = 0;
    domain_error("log_determinant_spd", "matrix argument", y,
                 "failed LDLT factorization");
  }

  // compute the inverse of A (needed for the derivative)
  m_d.setIdentity(m.rows(), m.cols());
  ldlt.solveInPlace(m_d);

  if (ldlt.isNegative() || (ldlt.vectorD().array() <= 1e-16).any()) {
    double y = 0;
    domain_error("log_determinant_spd", "matrix argument", y,
                 "matrix is negative definite");
  }

  double val = ldlt.vectorD().array().log().sum();

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
