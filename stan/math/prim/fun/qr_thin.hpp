#ifndef STAN_MATH_PRIM_FUN_QR_THIN_HPP
#define STAN_MATH_PRIM_FUN_QR_THIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>
#include <tuple>

namespace stan {
namespace math {

/**
 * Returns the thin QR decomposition
 *
 * @tparam EigMat type of the matrix
 * @param m Matrix.
 * @return A tuple containing (Q,R):
 * 1. Orthogonal matrix with minimal columns
 * 2. Upper triangular matrix with minimal rows
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
std::tuple<Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>,
           Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>>
qr_thin(const EigMat& m) {
  using matrix_t
      = Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>;
  if (unlikely(m.size() == 0)) {
    return std::make_tuple(matrix_t(0, 0), matrix_t(0, 0));
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
  matrix_t R = qr.matrixQR().topLeftCorner(min_size, m.cols());
  R.template triangularView<Eigen::StrictlyLower>().setZero();
  for (int i = 0; i < min_size; ++i) {
    if (R(i, i) < 0) {
      R.row(i) *= -1.0;
    }
  }
  return std::make_tuple(std::move(Q), std::move(R));
}

}  // namespace math
}  // namespace stan

#endif
