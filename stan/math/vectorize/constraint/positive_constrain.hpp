
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_POSITIVE_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_POSITIVE_CONSTRAIN_HPP
#include <stan/math/prim/constraint/positive_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return the positive value for the specified unconstrained input. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::EigenBase`, a `var_value` with inner type inheriting from
 * `Eigen::EigenBase`, a standard vector, or a scalar
 * @param x unconstrained value or container
 * @param[in, out] lp log density accumulator
 * @return positive constrained version of unconstrained value(s)
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto positive_constrain(const T& x, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      x, [&lp](auto&& v) { return positive_constrain<Jacobian>(v, lp); });
}


} // namespace math
} // namespace stan
#endif 

