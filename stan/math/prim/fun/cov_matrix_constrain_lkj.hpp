#ifndef STAN_MATH_PRIM_FUN_COV_MATRIX_CONSTRAIN_LKJ_HPP
#define STAN_MATH_PRIM_FUN_COV_MATRIX_CONSTRAIN_LKJ_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/corr_constrain.hpp>
#include <stan/math/prim/fun/positive_constrain.hpp>
#include <stan/math/prim/fun/read_cov_matrix.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the covariance matrix of the specified dimensionality
 * derived from constraining the specified vector of unconstrained
 * values.  The input vector must be of length \f$k \choose 2 +
 * k\f$.  The first \f$k \choose 2\f$ values in the input
 * represent unconstrained (partial) correlations and the last
 * \f$k\f$ are unconstrained standard deviations of the dimensions.
 *
 * <p>The transform scales the correlation matrix transform defined
 * in <code>corr_matrix_constrain(Matrix, size_t)</code>
 * with the constrained deviations.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations.
 * @param k Dimensionality of returned covariance matrix.
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cov_matrix_constrain_lkj(const T& x, size_t k) {
  size_t k_choose_2 = (k * (k - 1)) / 2;
  const auto& x_ref = to_ref(x);
  return read_cov_matrix(corr_constrain(x_ref.head(k_choose_2)),
                         positive_constrain(x_ref.tail(k)));
}

/**
 * Return the covariance matrix of the specified dimensionality
 * derived from constraining the specified vector of unconstrained
 * values and increment the specified log probability reference
 * with the log absolute Jacobian determinant.
 *
 * <p>The transform is defined as for
 * <code>cov_matrix_constrain(Matrix, size_t)</code>.
 *
 * <p>The log absolute Jacobian determinant is derived by
 * composing the log absolute Jacobian determinant for the
 * underlying correlation matrix as defined in
 * <code>cov_matrix_constrain(Matrix, size_t, T&)</code> with
 * the Jacobian of the transform of the correlation matrix
 * into a covariance matrix by scaling by standard deviations.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations.
 * @param k Dimensionality of returned covariance matrix.
 * @param lp Log probability reference to increment.
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cov_matrix_constrain_lkj(const T& x, size_t k, value_type_t<T>& lp) {
  size_t k_choose_2 = (k * (k - 1)) / 2;
  const auto& x_ref = x;
  return read_cov_matrix(corr_constrain(x_ref.head(k_choose_2)),
                         positive_constrain(x_ref.tail(k)), lp);
}

}  // namespace math
}  // namespace stan

#endif
