#ifndef STAN_MATH_PRIM_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/fmax.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform. This code is equivalent to
 * `softmax(sum_to_zero_constrain(y))`, but is more efficient and
 * stable if computed this way thanks to the use of the online
 * softmax algorithm courtesy of https://arxiv.org/abs/1805.02867.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Simplex of dimensionality K.
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr,
          require_not_st_var<Vec>* = nullptr>
inline plain_type_t<Vec> simplex_constrain(const Vec& y) {
  using T = value_type_t<Vec>;
  const auto N = y.size();

  plain_type_t<Vec> z = Eigen::VectorXd::Zero(N + 1);
  if (unlikely(N == 0)) {
    z.coeffRef(0) = 1;
    return z;
  }

  auto&& y_ref = to_ref(y);
  T sum_w(0);

  T d(0);  // sum of exponentials
  T max_val(0);
  T max_val_old(negative_infinity());

  for (int i = N; i > 0; --i) {
    double n = static_cast<double>(i);
    auto w = y_ref(i - 1) * inv_sqrt(n * (n + 1));
    sum_w += w;

    z.coeffRef(i - 1) += sum_w;
    z.coeffRef(i) -= w * n;

    max_val = fmax(max_val_old, z.coeff(i));
    d = d * exp(max_val_old - max_val) + exp(z.coeff(i) - max_val);
    max_val_old = max_val;
  }

  // above loop doesn't reach i==0
  max_val = fmax(max_val_old, z.coeff(0));
  d = d * exp(max_val_old - max_val) + exp(z.coeff(0) - max_val);

  z.array() = (z.array() - max_val).exp() / d;

  return z;
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform. This code is equivalent to
 * `softmax(sum_to_zero_constrain(y))`, but is more efficient and
 * stable if computed this way thanks to the use of the online
 * softmax algorithm courtesy of https://arxiv.org/abs/1805.02867.
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
inline plain_type_t<Vec> simplex_constrain(const Vec& y, Lp& lp) {
  using std::log;
  using T = value_type_t<Vec>;
  const auto N = y.size();

  plain_type_t<Vec> z = Eigen::VectorXd::Zero(N + 1);
  if (unlikely(N == 0)) {
    z.coeffRef(0) = 1;
    return z;
  }

  auto&& y_ref = to_ref(y);
  T sum_w(0);

  T d(0);  // sum of exponentials
  T max_val(0);
  T max_val_old(negative_infinity());

  for (int i = N; i > 0; --i) {
    double n = static_cast<double>(i);
    auto w = y_ref(i - 1) * inv_sqrt(n * (n + 1));
    sum_w += w;

    z.coeffRef(i - 1) += sum_w;
    z.coeffRef(i) -= w * n;

    max_val = fmax(max_val_old, z.coeff(i));
    d = d * exp(max_val_old - max_val) + exp(z.coeff(i) - max_val);
    max_val_old = max_val;
  }

  // above loop doesn't reach i==0
  max_val = fmax(max_val_old, z.coeff(0));
  d = d * exp(max_val_old - max_val) + exp(z.coeff(0) - max_val);

  z.array() = (z.array() - max_val).exp() / d;

  // equivalent to z.log().sum() + 0.5 * log(N + 1)
  lp += -(N + 1) * (max_val + log(d)) + 0.5 * log(N + 1);

  return z;
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
