
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#include <stan/math/prim/constraint/cholesky_corr_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return The cholesky of a `KxK` correlation matrix. If the `Jacobian`
 * parameter is `true`, the log density accumulator is incremented with the log
 * absolute Jacobian determinant of the transform.  All of the transforms are
 * specified with their Jacobians in the *Stan Reference Manual* chapter
 * Constraint Transforms.
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param y Linearly Serialized vector of size `(K * (K - 1))/2` holding the
 *  column major order elements of the lower triangurlar
 * @param K The size of the matrix to return
 * @param[in,out] lp log density accumulator
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cholesky_corr_constrain(const T& y, int K, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(y, [&lp, K](auto&& v) {
    return cholesky_corr_constrain<Jacobian>(v, K, lp);
  });
}


} // namespace math
} // namespace stan
#endif 

