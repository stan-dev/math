#ifndef STAN_MATH_PRIM_FUN_QR_HPP
#define STAN_MATH_PRIM_FUN_QR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>
#include <tuple>

namespace stan {
namespace math {

/**
 * Returns the fat QR decomposition
 *
 * @tparam EigMat type of the matrix
 * @param m Matrix.
 * @return A tuple containing (Q,R):
 * 1. Orthogonal matrix with maximal columns
 * 2. Upper triangular matrix with maximal rows
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
std::tuple<Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>,
           Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>>
qr(const EigMat& m) {
  using matrix_t
      = Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>;
  check_nonzero_size("qr", "m", m);
  Eigen::HouseholderQR<matrix_t> qr(m.rows(), m.cols());
  qr.compute(m);
  matrix_t Q = qr.householderQ();
  const int min_size = std::min(m.rows(), m.cols());
  for (int i = 0; i < min_size; i++) {
    if (qr.matrixQR().coeff(i, i) < 0) {
      Q.col(i) *= -1.0;
    }
  }
  matrix_t R = qr.matrixQR();
  if (m.rows() > m.cols()) {
    R.bottomRows(m.rows() - m.cols()).setZero();
  }
  for (int i = 0; i < min_size; i++) {
    for (int j = 0; j < i; j++) {
      R.coeffRef(i, j) = 0.0;
    }
    if (R(i, i) < 0) {
      R.row(i) *= -1.0;
    }
  }
  return std::make_tuple(std::move(Q), std::move(R));
}

}  // namespace math
}  // namespace stan

#endif
