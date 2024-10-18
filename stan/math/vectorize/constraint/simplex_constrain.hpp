
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP
#include <stan/math/prim/constraint/simplex_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return the simplex corresponding to the specified free vector. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam Vec A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param[in] y free vector
 * @param[in, out] lp log density accumulator
 * @return simplex of dimensionality one greater than `y`
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto simplex_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      y, [&lp](auto&& v) { return simplex_constrain<Jacobian>(v, lp); });
}


} // namespace math
} // namespace stan
#endif 

