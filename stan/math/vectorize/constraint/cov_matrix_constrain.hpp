
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_COV_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_COV_MATRIX_CONSTRAIN_HPP
#include <stan/math/prim/constraint/cov_matrix_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return the symmetric, positive-definite matrix of dimensions K by K resulting
 * from transforming the specified finite vector of size K plus (K choose 2). If
 * the `Jacobian` parameter is `true`, the log density accumulator is
 * incremented with the log absolute Jacobian determinant of the transform.  All
 * of the transforms are specified with their Jacobians in the *Stan Reference
 * Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param x The vector to convert to a covariance matrix
 * @param K The dimensions of the resulting covariance matrix
 * @param[in, out] lp log density accumulator
 * @throws std::domain_error if (x.size() != K + (K choose 2)).
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cov_matrix_constrain(const T& x, Eigen::Index K,
                                 return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(x, [&lp, K](auto&& v) {
    return cov_matrix_constrain<Jacobian>(v, K, lp);
  });
}


} // namespace math
} // namespace stan
#endif 

