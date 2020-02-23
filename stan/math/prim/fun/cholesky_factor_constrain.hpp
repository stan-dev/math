#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector.  A total of (N choose 2) + N + (M - N) * N
 * elements are required to read an M by N Cholesky factor.
 *
 * @tparam T type of elements in the matrix
 * @param x Vector of unconstrained values
 * @param M Number of rows
 * @param N Number of columns
 * @return Cholesky factor
 */
template <typename EigVec, typename = require_eigen_t<EigVec>>
auto cholesky_factor_constrain(EigVec&& x, int M, int N) {
  using eigen_scalar = value_type_t<EigVec>;
  using std::exp;
  check_greater_or_equal("cholesky_factor_constrain",
                         "num rows (must be greater or equal to num cols)", M,
                         N);
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  Eigen::Matrix<eigen_scalar, Eigen::Dynamic, Eigen::Dynamic> y(M, N);
  int kk = 0;
  // NOTE: Why does x come in row major order?
  for (auto i = 0; i < N; i++) {
    y.row(i).segment(0, i + 1) = x.segment(i + kk, i + 1);
    kk += i;
  }
  y.template triangularView<Eigen::StrictlyUpper>().setZero();
  y.diagonal().array() = y.diagonal().array().exp();
  return y;
}

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector and increment the specified log probability
 * reference with the log Jacobian adjustment of the transform.  A total
 * of (N choose 2) + N + N * (M - N) free parameters are required to read
 * an M by N Cholesky factor.
 *
 * @tparam T type of elements in the matrix
 * @param x Vector of unconstrained values
 * @param M Number of rows
 * @param N Number of columns
 * @param lp Log probability that is incremented with the log Jacobian
 * @return Cholesky factor
 */
template <typename EigMat, typename Index1, typename Index2, typename T, typename = require_eigen_t<EigMat>>
auto cholesky_factor_constrain(const EigMat& xx, Index1 M, Index2 N, T& lp) {
  using eigen_scalar = value_type_t<EigMat>;
  const Eigen::Ref<const typename EigMat::PlainObject>& x = xx;
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  int pos = 0;
  Eigen::Matrix<eigen_scalar, -1, 1> log_jacobians(N);
  for (auto n = 0; n < N; ++n) {
    pos += n;
    log_jacobians[n] = x.coeffRef(pos);
    pos++;
  }
  lp += sum(log_jacobians);
  return cholesky_factor_constrain(std::forward<EigMat>(x), M, N);
}

}  // namespace math
}  // namespace stan

#endif
