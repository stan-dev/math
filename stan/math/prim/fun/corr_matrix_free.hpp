#ifndef STAN_MATH_PRIM_FUN_CORR_MATRIX_FREE_HPP
#define STAN_MATH_PRIM_FUN_CORR_MATRIX_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/factor_cov_matrix.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the vector of unconstrained partial correlations that
 * define the specified correlation matrix when transformed.
 *
 * <p>The constraining transform is defined as for
 * <code>corr_matrix_constrain(Matrix, size_t)</code>.  The
 * inverse transform in this function is simpler in that it only
 * needs to compute the \f$k \choose 2\f$ partial correlations
 * and then free those.
 *
 * @tparam T type of scalar
 * @param y The correlation matrix to free.
 * @return Vector of unconstrained values that produce the
 * specified correlation matrix when transformed.
 * @throw std::domain_error if the correlation matrix has no
 *    elements or is not a square matrix.
 * @throw std::runtime_error if the correlation matrix cannot be
 *    factorized by factor_cov_matrix() or if the sds returned by
 *    factor_cov_matrix() on log scale are unconstrained.
 */
template <typename EigMat, typename = require_eigen_t<EigMat>>
inline auto corr_matrix_free(EigMat&& y) {
  check_square("corr_matrix_free", "y", y);
  check_nonzero_size("corr_matrix_free", "y", y);

  using Eigen::Array;
  using Eigen::Dynamic;
  using eigen_scalar = value_type_t<EigMat>;

  auto k = y.rows();
  auto k_choose_2 = (k * (k - 1)) / 2;
  auto ret_vals = factor_cov_matrix(std::forward<EigMat>(y));
  check_bounded("corr_matrix_free", "log(sd)", std::get<1>(ret_vals),
   -CONSTRAINT_TOLERANCE, CONSTRAINT_TOLERANCE);
  // NOTE: removing eval here causes vector of zeros to return (???)
  return std::get<0>(ret_vals).matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
