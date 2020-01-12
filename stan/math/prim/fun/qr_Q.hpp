#ifndef STAN_MATH_PRIM_FUN_QR_Q_HPP
#define STAN_MATH_PRIM_FUN_QR_Q_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Returns the orthogonal factor of the fat QR decomposition
 *
 * @tparam T type of elements in the matrix
 * @param m Matrix.
 * @return Orthogonal matrix with maximal columns
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> qr_Q(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  check_nonzero_size("qr_Q", "m", m);
  Eigen::HouseholderQR<matrix_t> qr(m.rows(), m.cols());
  qr.compute(m);
  matrix_t Q = qr.householderQ();
  const int min_size = std::min(m.rows(), m.cols());
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
