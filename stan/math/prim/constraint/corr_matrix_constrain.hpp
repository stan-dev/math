#ifndef STAN_MATH_PRIM_CONSTRAINT_CORR_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CORR_MATRIX_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/corr_constrain.hpp>
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
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x Vector of unconstrained partial correlations.
 * @param k Dimensionality of returned correlation matrix.
 * @throw std::invalid_argument if x is not a valid correlation
 * matrix.
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
corr_matrix_constrain(const T& x, Eigen::Index k) {
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
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x Vector of unconstrained partial correlations.
 * @param k Dimensionality of returned correlation matrix.
 * @param lp Log probability reference to increment.
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
corr_matrix_constrain(const T& x, Eigen::Index k, return_type_t<T>& lp) {
  Eigen::Index k_choose_2 = (k * (k - 1)) / 2;
  check_size_match("cov_matrix_constrain", "x.size()", x.size(), "k_choose_2",
                   k_choose_2);
  return read_corr_matrix(corr_constrain(x, lp), k, lp);
}

/**
 * Return the correlation matrix of the specified dimensionality derived from
 * the specified vector of unconstrained values. The input vector must be of
 * length \f${k \choose 2} = \frac{k(k-1)}{2}\f$.  The values in the input
 * vector represent unconstrained (partial) correlations among the dimensions.
 * If the `Jacobian` parameter is `true`, the log density accumulator is
 * incremented with the log absolute Jacobian determinant of the transform.  All
 * of the transforms are specified with their Jacobians in the *Stan Reference
 * Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @param x Vector of unconstrained partial correlations
 * @param k Dimensionality of returned correlation matrix
 * @param[in,out] lp log density accumulator
 */
template <bool Jacobian, typename T, require_not_std_vector_t<T>* = nullptr>
inline auto corr_matrix_constrain(const T& x, Eigen::Index k,
                                  return_type_t<T>& lp) {
  if (Jacobian) {
    return corr_matrix_constrain(x, k, lp);
  } else {
    return corr_matrix_constrain(x, k);
  }
}

/**
 * Return the correlation matrix of the specified dimensionality derived from
 * the specified vector of unconstrained values. The input vector must be of
 * length \f${k \choose 2} = \frac{k(k-1)}{2}\f$.  The values in the input
 * vector represent unconstrained (partial) correlations among the dimensions.
 * If the `Jacobian` parameter is `true`, the log density accumulator is
 * incremented with the log absolute Jacobian determinant of the transform.  All
 * of the transforms are specified with their Jacobians in the *Stan Reference
 * Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param y Vector of unconstrained partial correlations
 * @param K Dimensionality of returned correlation matrix
 * @param[in,out] lp log density accumulator
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto corr_matrix_constrain(const T& y, int K, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(y, [&lp, K](auto&& v) {
    return corr_matrix_constrain<Jacobian>(v, K, lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
