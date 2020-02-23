#ifndef STAN_MATH_PRIM_FUN_COV_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_COV_MATRIX_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the symmetric, positive-definite matrix of dimensions K
 * by K resulting from transforming the specified finite vector of
 * size K plus (K choose 2).
 *
 * <p>See <code>cov_matrix_free()</code> for the inverse transform.
 *
 * @tparam T type of elements in the vector
 * @param x The vector to convert to a covariance matrix.
 * @param K The number of rows and columns of the resulting
 * covariance matrix.
 * @throws std::invalid_argument if (x.size() != K + (K choose 2)).
 */
template <typename EigMat, typename Index, typename = require_eigen_t<EigMat>>
auto cov_matrix_constrain(EigMat&& x, Index K) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using eigen_scalar = value_type_t<EigMat>;
  check_size_match("cov_matrix_constrain", "x.size()", x.size(),
                   "K + (K choose 2)", (K * (K + 1)) / 2);
  Eigen::Matrix<eigen_scalar, Dynamic, Dynamic> L(K, K);
  int kk = 0;
  // NOTE: Why does x come in row major order?
  for (Index i = 0; i < K; i++) {
    L.row(i).segment(0, i + 1) = x.segment(i + kk, i + 1);
    kk += i;
  }
  L.diagonal().array() = L.diagonal().array().exp();
  L.template triangularView<Eigen::StrictlyUpper>().setZero();
  return multiply_lower_tri_self_transpose(L);
}

/**
 * Return the symmetric, positive-definite matrix of dimensions K
 * by K resulting from transforming the specified finite vector of
 * size K plus (K choose 2).
 *
 * <p>See <code>cov_matrix_free()</code> for the inverse transform.
 *
 * @tparam T type of elements in the vector
 * @param x The vector to convert to a covariance matrix.
 * @param K The dimensions of the resulting covariance matrix.
 * @param lp Reference
 * @throws std::domain_error if (x.size() != K + (K choose 2)).
 */
template <typename EigMat, typename Index, typename T,
          typename = require_eigen_t<EigMat>>
auto cov_matrix_constrain(EigMat&& x, Index K, T& lp) {
  // TODO(Restrict these?)
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using std::log;
  using eigen_scalar = value_type_t<EigMat>;
  check_size_match("cov_matrix_constrain", "x.size()", x.size(),
                   "K + (K choose 2)", (K * (K + 1)) / 2);
  Eigen::Matrix<eigen_scalar, Dynamic, Dynamic> L(K, K);
  int kk = 0;
  // NOTE: x comes in with row major ordered lower triangular as vector
  for (Index i = 0; i < K; i++) {
    L.row(i).segment(0, i + 1) = x.segment(i + kk, i + 1);
    kk += i;
  }
  L.diagonal().array() = L.diagonal().array().exp();
  L.template triangularView<Eigen::StrictlyUpper>().setZero();
  // Jacobian for complete transform, including exp() above
  lp += (K * LOG_TWO);  // needless constant; want propto
  lp += ((K - Eigen::ArrayXd::LinSpaced(K, 1, K)) * log(L.diagonal()).array())
            .sum();  // only +1 because index from 0
  return multiply_lower_tri_self_transpose(L);
}

}  // namespace math
}  // namespace stan

#endif
