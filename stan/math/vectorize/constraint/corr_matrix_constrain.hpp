
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_CORR_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_CORR_MATRIX_CONSTRAIN_HPP
#include <stan/math/prim/constraint/corr_matrix_constrain.hpp>
namespace stan {
namespace math {

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


} // namespace math
} // namespace stan
#endif 

