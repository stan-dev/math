#ifndef STAN_MATH_PRIM_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The transform is based on a centered stick-breaking process.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Simplex of dimensionality K.
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr,
          require_not_st_var<Vec>* = nullptr>
inline plain_type_t<Vec> simplex_constrain(Vec&& y) {
  // cut & paste simplex_constrain(Eigen::Matrix, T) w/o Jacobian
  using std::log;
  using T = value_type_t<Vec>;

  auto&& y_ref = to_ref(std::forward<Vec>(y));
  const Eigen::Index Km1 = y_ref.size();
  plain_type_t<Vec> x(Km1 + 1);
  T stick_len(1.0);
  for (Eigen::Index k = 0; k < Km1; ++k) {
    T z_k = inv_logit(y_ref.coeff(k) - log(Km1 - k));
    x.coeffRef(k) = stick_len * z_k;
    stick_len -= x.coeff(k);
  }
  x.coeffRef(Km1) = stick_len;
  return x;
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam Vec type of the vector
 * @tparam Lp A scalar type for the lp argument. The scalar type of Vec should
 * be convertable to this.
 * @param y Free vector input of dimensionality K - 1.
 * @param lp Log probability reference to increment.
 * @return Simplex of dimensionality K.
 */
template <typename Vec, typename Lp, require_eigen_vector_t<Vec>* = nullptr,
          require_not_st_var<Vec>* = nullptr,
          require_convertible_t<value_type_t<Vec>, Lp>* = nullptr>
inline plain_type_t<Vec> simplex_constrain(Vec&& y, Lp& lp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::log;
  using T = value_type_t<Vec>;
  auto&& y_ref = to_ref(std::forward<Vec>(y));
  const Eigen::Index Km1 = y.size();  // K = Km1 + 1
  plain_type_t<Vec> x(Km1 + 1);
  T stick_len(1.0);
  for (Eigen::Index k = 0; k < Km1; ++k) {
    double eq_share = -log(Km1 - k);  // = logit(1.0/(Km1 + 1 - k));
    T adj_y_k = y_ref.coeff(k) + eq_share;
    T z_k = inv_logit(adj_y_k);
    x.coeffRef(k) = stick_len * z_k;
    lp += log(stick_len) - log1p_exp(-adj_y_k) - log1p_exp(adj_y_k);
    stick_len -= x.coeff(k);  // equivalently *= (1 - z_k);
  }
  x.coeffRef(Km1) = stick_len;  // no Jacobian contrib for last dim
  return x;
}

/**
 * Return the simplex corresponding to the specified free vector.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam Vec A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param[in] y free vector
 * @return simplex of dimensionality one greater than `y`
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto simplex_constrain(T&& y) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [](auto&& v) {
    return simplex_constrain(std::forward<decltype(v)>(v));
  });
}

/**
 * Return the simplex corresponding to the specified free vector.
 * This overload handles looping over the elements of a standard vector.
 *
 * @tparam Vec A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @tparam Lp Scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param[in] y free vector
 * @param[in, out] lp log density accumulator
 * @return simplex of dimensionality one greater than `y`
 */
template <typename T, typename Lp, require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto simplex_constrain(T&& y, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [&lp](auto&& v) {
    return simplex_constrain(std::forward<decltype(v)>(v), lp);
  });
}

/**
 * Return the simplex corresponding to the specified free vector. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam Vec A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @tparam Lp A scalar type for the lp argument. The scalar type of Vec should
 * be convertable to this.
 * @param[in] y free vector
 * @param[in, out] lp log density accumulator
 * @return simplex of dimensionality one greater than `y`
 */
template <bool Jacobian, typename Vec, typename Lp,
          require_convertible_t<return_type_t<Vec>, Lp>* = nullptr>
inline plain_type_t<Vec> simplex_constrain(Vec&& y, Lp& lp) {
  if constexpr (Jacobian) {
    return simplex_constrain(std::forward<Vec>(y), lp);
  } else {
    return simplex_constrain(std::forward<Vec>(y));
  }
}

}  // namespace math
}  // namespace stan

#endif
