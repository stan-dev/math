
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_COV_MATRIX_CONSTRAIN_LKJ_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_COV_MATRIX_CONSTRAIN_LKJ_HPP
#include <stan/math/prim/constraint/cov_matrix_constrain_lkj.hpp>
namespace stan {
namespace math {

/**
 * Return the covariance matrix of the specified dimensionality derived from
 * constraining the specified vector of unconstrained values. If the `Jacobian`
 * parameter is `true`, the log density accumulator is incremented with the log
 * absolute Jacobian determinant of the transform.  All of the transforms are
 * specified with their Jacobians in the *Stan Reference Manual* chapter
 * Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time rows or columns equal to 1
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations
 * @param k Dimensionality of returned covariance matrix
 * @param[in, out] lp log density accumulator
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cov_matrix_constrain_lkj(const T& x, size_t k,
                                     return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(x, [&lp, k](auto&& v) {
    return cov_matrix_constrain_lkj<Jacobian>(v, k, lp);
  });
}


} // namespace math
} // namespace stan
#endif 

