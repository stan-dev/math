
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_UNIT_VECTOR_CONSTRAIN_HPP
#include <stan/math/prim/constraint/unit_vector_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector y. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param y vector of K unrestricted variables
 * @param[in, out] lp log density accumulator
 * @return Unit length vector of dimension K
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto unit_vector_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      y, [&lp](auto&& v) { return unit_vector_constrain<Jacobian>(v, lp); });
}


} // namespace math
} // namespace stan
#endif 

