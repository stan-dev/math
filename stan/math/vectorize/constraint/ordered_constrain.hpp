
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_ORDERED_CONSTRAIN_HPP
#include <stan/math/prim/constraint/ordered_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return a positive valued, increasing ordered vector derived from the
 * specified free vector. The returned constrained vector will have the same
 * dimensionality as the specified free vector. If the `Jacobian` parameter is
 * `true`, the log density accumulator is incremented with the log absolute
 * Jacobian determinant of the transform. All of the transforms are specified
 * with their Jacobians in the *Stan Reference Manual* chapter Constraint
 * Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param x Free vector of scalars
 * @param[in, out] lp log density accumulator
 * @return Positive, increasing ordered vector.
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto ordered_constrain(const T& x, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      x, [&lp](auto&& v) { return ordered_constrain<Jacobian>(v, lp); });
}


} // namespace math
} // namespace stan
#endif 

