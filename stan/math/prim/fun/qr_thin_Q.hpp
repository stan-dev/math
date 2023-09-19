#ifndef STAN_MATH_PRIM_FUN_QR_THIN_Q_HPP
#define STAN_MATH_PRIM_FUN_QR_THIN_Q_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns the orthogonal factor of the thin QR decomposition
 *
 * @tparam EigMat type of the matrix
 * @param m Matrix.
 * @return Orthogonal matrix with minimal columns
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic> qr_thin_Q(
    const EigMat& m) {
  using matrix_t
      = Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>;
  if (unlikely(m.size() == 0)) {
    return matrix_t(0, 0);
  }
  Eigen::HouseholderQR<matrix_t> qr(m.rows(), m.cols());
  qr.compute(m);
  const int min_size = std::min(m.rows(), m.cols());
  matrix_t Q = qr.householderQ() * matrix_t::Identity(m.rows(), min_size);
  for (int i = 0; i < min_size; i++) {
    if (qr.matrixQR().coeff(i, i) < 0) {
      Q.col(i) *= -1.0;
    }
  }
  return Q;
}

}  // namespace math
}  // namespace stan

#endif
