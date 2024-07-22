#ifndef STAN_MATH_PRIM_CONSTRAINT_POSITIVE_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_POSITIVE_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an increasing positive ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam T type of elements in the vector
 * @param x Free vector of scalars.
 * @return Positive, increasing ordered vector.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr,
          require_not_st_var<EigVec>* = nullptr>
inline auto positive_ordered_constrain(const EigVec& x) {
  using std::exp;
  Eigen::Index k = x.size();
  plain_type_t<EigVec> y(k);
  if (k == 0) {
    return y;
  }
  const auto& x_ref = to_ref(x);
  y.coeffRef(0) = exp(x_ref.coeff(0));
  for (Eigen::Index i = 1; i < k; ++i) {
    y.coeffRef(i) = y.coeff(i - 1) + exp(x_ref.coeff(i));
  }
  return y;
}

/**
 * Return a positive valued, increasing positive ordered vector derived
 * from the specified free vector and increment the specified log
 * probability reference with the log absolute Jacobian determinant
 * of the transform.  The returned constrained vector
 * will have the same dimensionality as the specified free vector.
 *
 * @tparam T type of elements in the vector
 * @param x Free vector of scalars.
 * @param lp Log probability reference.
 * @return Positive, increasing ordered vector.
 */
template <typename Vec, require_col_vector_t<Vec>* = nullptr>
inline auto positive_ordered_constrain(const Vec& x, return_type_t<Vec>& lp) {
  const auto& x_ref = to_ref(x);
  lp += sum(x_ref);
  return positive_ordered_constrain(x_ref);
}

/**
 * Return a positive valued, increasing positive ordered vector derived from the
 * specified free vector. The returned constrained vector will have the same
 * dimensionality as the specified free vector. If the `Jacobian` parameter is
 * `true`, the log density accumulator is incremented with the log absolute
 * Jacobian determinant of the transform.  All of the transforms are specified
 * with their Jacobians in the *Stan Reference Manual* chapter Constraint
 * Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam Vec A type inheriting from `Eigen::EigenBase`, a `var_value` with
 * inner type inheriting from `Eigen::EigenBase`
 * @param x Free vector of scalars
 * @param[in, out] lp log density accumulato
 * @return Positive, increasing ordered vector
 */
template <bool Jacobian, typename Vec, require_not_std_vector_t<Vec>* = nullptr>
inline auto positive_ordered_constrain(const Vec& x, return_type_t<Vec>& lp) {
  if (Jacobian) {
    return positive_ordered_constrain(x, lp);
  } else {
    return positive_ordered_constrain(x);
  }
}

/**
 * Return a positive valued, increasing positive ordered vector derived from the
 * specified free vector. The returned constrained vector will have the same
 * dimensionality as the specified free vector. If the `Jacobian` parameter is
 * `true`, the log density accumulator is incremented with the log absolute
 * Jacobian determinant of the transform.  All of the transforms are specified
 * with their Jacobians in the *Stan Reference Manual* chapter Constraint
 * Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam Vec A standard vector with inner type inheriting from
 * `Eigen::EigenBase`, a `var_value` with inner type inheriting from
 * `Eigen::EigenBase`
 * @param x Free vector of scalars
 * @param[in, out] lp log density accumulato
 * @return Positive, increasing ordered vector
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto positive_ordered_constrain(const T& x, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(x, [&lp](auto&& v) {
    return positive_ordered_constrain<Jacobian>(v, lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
