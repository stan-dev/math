#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector y.
 *
 * See <a
 * href="https://en.wikipedia.org/wiki/N-sphere#Generating_random_points">the
 * Wikipedia page on generating random points on an N-sphere</a>.
 *
 * @tparam T type inheriting from `EigenBase` that does not have an fvar
 *  scalar type.
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr,
          require_not_vt_autodiff<T>* = nullptr>
inline plain_type_t<T> unit_vector_constrain(const T& y) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "y", y);
  auto&& y_ref = to_ref(y);
  value_type_t<T> SN = dot_self(y_ref);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  return y_ref.array() / sqrt(SN);
}

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam T1 type inheriting from `EigenBase` that does not have an fvar
 *  scalar type.
 *
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 * @param lp Log probability reference to increment.
 */
template <typename T1, typename T2, require_eigen_col_vector_t<T1>* = nullptr,
          require_all_not_vt_autodiff<T1, T2>* = nullptr>
inline plain_type_t<T1> unit_vector_constrain(const T1& y, T2& lp) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "y", y);
  auto&& y_ref = to_ref(y);
  value_type_t<T1> SN = dot_self(y_ref);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  lp -= 0.5 * SN;
  return y_ref.array() / sqrt(SN);
}

/**
 * Return the unit length vector corresponding to the free vector y. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @param y vector of K unrestricted variables
 * @param[in, out] lp log density accumulator
 * @return Unit length vector of dimension K
 */
template <bool Jacobian, typename T, require_not_std_vector_t<T>* = nullptr>
inline auto unit_vector_constrain(const T& y, return_type_t<T>& lp) {
  if (Jacobian) {
    return unit_vector_constrain(y, lp);
  } else {
    return unit_vector_constrain(y);
  }
}

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

}  // namespace math
}  // namespace stan

#endif
