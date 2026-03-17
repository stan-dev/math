#ifndef STAN_MATH_PRIM_CONSTRAINT_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an increasing ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam T type of the vector
 * @param x Free vector of scalars.
 * @return Positive, increasing ordered vector.
 * @tparam T Type of scalar.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr,
          require_not_st_var<EigVec>* = nullptr>
inline plain_type_t<EigVec> ordered_constrain(EigVec&& x) {
  using std::exp;
  auto&& x_ref = to_ref(std::forward<EigVec>(x));
  Eigen::Index k = x_ref.size();
  plain_type_t<EigVec> y(k);
  if (unlikely(k == 0)) {
    return y;
  }
  y[0] = x_ref[0];
  for (Eigen::Index i = 1; i < k; ++i) {
    y.coeffRef(i) = y.coeff(i - 1) + exp(x_ref.coeff(i));
  }
  return y;
}

/**
 * Return a positive valued, increasing ordered vector derived
 * from the specified free vector and increment the specified log
 * probability reference with the log absolute Jacobian determinant
 * of the transform.  The returned constrained vector
 * will have the same dimensionality as the specified free vector.
 *
 * @tparam T type of the vector
 * @tparam Lp A scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param x Free vector of scalars.
 * @param lp Log probability reference.
 * @return Positive, increasing ordered vector.
 */
template <typename EigVec, typename Lp,
          require_eigen_col_vector_t<EigVec>* = nullptr,
          require_convertible_t<value_type_t<EigVec>, Lp>* = nullptr>
inline auto ordered_constrain(EigVec&& x, Lp& lp) {
  auto&& x_ref = to_ref(std::forward<EigVec>(x));
  if (likely(x_ref.size() > 1)) {
    lp += sum(x_ref.tail(x_ref.size() - 1));
  }
  return ordered_constrain(std::forward<decltype(x_ref)>(x_ref));
}

/**
 * Return a positive valued, increasing ordered vector derived from the
 * specified free vector. The returned constrained vector will have the same
 * dimensionality as the specified free vector.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param x Free vector of scalars
 * @return Positive, increasing ordered vector.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto ordered_constrain(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return ordered_constrain(std::forward<decltype(v)>(v));
  });
}

/**
 * Return a positive valued, increasing ordered vector derived from the
 * specified free vector. The returned constrained vector will have the same
 * dimensionality as the specified free vector.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @tparam Lp Scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param x Free vector of scalars
 * @param[in, out] lp log density accumulator or empty
 * @return Positive, increasing ordered vector.
 */
template <typename T, typename Lp, require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto ordered_constrain(T&& x, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [&lp](auto&& v) {
    return ordered_constrain(std::forward<decltype(v)>(v), lp);
  });
}

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
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column, or a standard vector thereof
 * @tparam Lp A scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param x Free vector of scalars
 * @param[in, out] lp log density accumulator
 * @return Positive, increasing ordered vector.
 */
template <bool Jacobian, typename T, typename Lp,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto ordered_constrain(T&& x, Lp& lp) {
  if constexpr (Jacobian) {
    return ordered_constrain(std::forward<T>(x), lp);
  } else {
    return ordered_constrain(std::forward<T>(x));
  }
}

}  // namespace math
}  // namespace stan

#endif
