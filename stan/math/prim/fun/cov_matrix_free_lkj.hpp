#ifndef STAN_MATH_PRIM_FUN_COV_MATRIX_FREE_LKJ_HPP
#define STAN_MATH_PRIM_FUN_COV_MATRIX_FREE_LKJ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the vector of unconstrained partial correlations and
 * deviations that transform to the specified covariance matrix.
 *
 * <p>The constraining transform is defined as for
 * <code>cov_matrix_constrain(Matrix, size_t)</code>.  The
 * inverse first factors out the deviations, then applies the
 * freeing transform of <code>corr_matrix_free(Matrix&)</code>.
 *
 * @tparam T type of elements in the matrix
 * @param y Covariance matrix to free.
 * @return Vector of unconstrained values that transforms to the
 * specified covariance matrix.
 * @throw std::domain_error if the correlation matrix has no
 *    elements or is not a square matrix.
 * @throw std::runtime_error if the correlation matrix cannot be
 *    factorized by factor_cov_matrix()
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> cov_matrix_free_lkj(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& y) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using size_type = typename index_type<Matrix<T, Dynamic, Dynamic>>::type;

  check_nonzero_size("cov_matrix_free_lkj", "y", y);
  check_square("cov_matrix_free_lkj", "y", y);
  size_type k = y.rows();
  size_type k_choose_2 = (k * (k - 1)) / 2;
  Array<T, Dynamic, 1> cpcs(k_choose_2);
  Array<T, Dynamic, 1> sds(k);
  bool successful = factor_cov_matrix(y, cpcs, sds);
  if (!successful) {
    throw_domain_error("cov_matrix_free_lkj", "factor_cov_matrix failed on y",
                       "", "");
  }
  Matrix<T, Dynamic, 1> x(k_choose_2 + k);
  size_type pos = 0;
  for (size_type i = 0; i < k_choose_2; ++i) {
    x[pos++] = cpcs[i];
  }
  for (size_type i = 0; i < k; ++i) {
    x[pos++] = sds[i];
  }
  return x;
}

}  // namespace math
}  // namespace stan

#endif
