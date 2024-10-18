
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_STOCHASTIC_COLUMN_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_STOCHASTIC_COLUMN_CONSTRAIN_HPP
#include <stan/math/prim/constraint/stochastic_column_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return a vector of column stochastic matrices. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic columns
 * @param[in] y free vector
 * @param[in, out] lp log density accumulator
 * @return Standard vector containing matrices with simplex columns of
 * dimensionality (K, M).
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_column_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(y, [&lp](auto&& v) {
    return stochastic_column_constrain<Jacobian>(v, lp);
  });
}


} // namespace math
} // namespace stan
#endif 

