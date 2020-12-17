#ifndef STAN_MATH_REV_FUN_CORR_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_CORR_MATRIX_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/corr_constrain.hpp>
#include <stan/math/prim/fun/read_corr_matrix.hpp>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Return the correlation matrix of the specified dimensionality
 * derived from the specified vector of unconstrained values.  The
 * input vector must be of length \f${k \choose 2} =
 * \frac{k(k-1)}{2}\f$.  The values in the input vector represent
 * unconstrained (partial) correlations among the dimensions.
 *
 * <p>The transform based on partial correlations is as specified
 * in
 *
 * <ul><li> Lewandowski, Daniel, Dorota Kurowicka, and Harry
 * Joe. 2009.  Generating random correlation matrices based on
 * vines and extended onion method.  <i>Journal of Multivariate
 * Analysis</i> <b>100</b>:1989–-2001.  </li></ul>
 *
 * <p>The free vector entries are first constrained to be
 * valid correlation values using <code>corr_constrain(T)</code>.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x Vector of unconstrained partial correlations.
 * @param k Dimensionality of returned correlation matrix.
 * @throw std::invalid_argument if x is not a valid correlation
 * matrix.
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> corr_matrix_constrain(const T& x, Eigen::Index k) {
  Eigen::Index k_choose_2 = (k * (k - 1)) / 2;
  check_size_match("cov_matrix_constrain", "x.size()", x.size(), "k_choose_2",
                   k_choose_2);
  return read_corr_matrix(corr_constrain(x), k);
}

/**
 * Return the correlation matrix of the specified dimensionality
 * derived from the specified vector of unconstrained values.  The
 * input vector must be of length \f${k \choose 2} =
 * \frac{k(k-1)}{2}\f$.  The values in the input vector represent
 * unconstrained (partial) correlations among the dimensions.
 *
 * <p>The transform is as specified for
 * <code>corr_matrix_constrain(Matrix, size_t)</code>; the
 * paper it cites also defines the Jacobians for correlation inputs,
 * which are composed with the correlation constrained Jacobians
 * defined in <code>corr_constrain(T, double)</code> for
 * this function.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x Vector of unconstrained partial correlations.
 * @param k Dimensionality of returned correlation matrix.
 * @param lp Log probability reference to increment.
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> corr_matrix_constrain(const T& x, Eigen::Index k,
                                                 scalar_type_t<T>& lp) {
  Eigen::Index k_choose_2 = (k * (k - 1)) / 2;
  check_size_match("cov_matrix_constrain", "x.size()", x.size(), "k_choose_2",
                   k_choose_2);
  return read_corr_matrix(corr_constrain(x, lp), k, lp);
}

}  // namespace math
}  // namespace stan

#endif
