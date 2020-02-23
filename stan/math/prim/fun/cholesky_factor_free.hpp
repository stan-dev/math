#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_FREE_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_FREE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Return the unconstrained vector of parameters correspdonding to
 * the specified Cholesky factor.  A Cholesky factor must be lower
 * triangular and have positive diagonal elements.
 *
 * @tparam T type of elements in the Cholesky factor
 * @param y Cholesky factor.
 * @return Unconstrained parameters for Cholesky factor.
 * @throw std::domain_error If the matrix is not a Cholesky factor.
 */
template <typename EigMat, typename = require_eigen_t<EigMat>>
auto cholesky_factor_free(EigMat&& y) {
  using eigen_scalar = value_type_t<EigMat>;
  using std::log;
  check_cholesky_factor("cholesky_factor_free", "y", y);
  int M = y.rows();
  int N = y.cols();
  Eigen::Matrix<eigen_scalar, Eigen::Dynamic, 1> x((N * (N + 1)) / 2 + (M - N) * N);
  // NOTE: Couldn't figure out indexing for col major :-(
  int kk = 0;
  for (auto i = 0; i < M; i++) {
    x.segment(i + kk, i + 1) = y.row(i).segment(0, i + 1);
    kk += i;
    x(i + kk) = log(x(i + kk));
  }
  return x;
}

}  // namespace math
}  // namespace stan

#endif
