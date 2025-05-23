#ifndef STAN_MATH_PRIM_CONSTRAINT_OFFSET_MULTIPLIER_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_OFFSET_MULTIPLIER_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/constraint/identity_constrain.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * <p>If the offset is zero and the multiplier is one this
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T A scalar type or type inheriting from `Eigen::DenseBase`
 * @tparam M A scalar type or type inheriting from `Eigen::DenseBase`
 * @tparam S A scalar type or type inheriting from `Eigen::DenseBase`
 * @param[in] x Unconstrained scalar input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T, M, S>* = nullptr,
          require_all_not_std_vector_t<T, M, S>* = nullptr>
inline auto offset_multiplier_constrain(T&& x, M&& mu, S&& sigma) {
  if (is_matrix_v<T> && is_matrix<M>::value) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "mu", mu);
  }
  if (is_matrix_v<T> && is_matrix_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "sigma", sigma);
  } else if (is_matrix<M>::value && is_matrix_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "mu", mu, "sigma",
                        sigma);
  }
  auto&& mu_ref = to_ref(std::forward<M>(mu));
  auto&& sigma_ref = to_ref(std::forward<S>(sigma));
  check_finite("offset_multiplier_constrain", "offset", value_of_rec(mu_ref));
  check_positive_finite("offset_multiplier_constrain", "multiplier",
                        value_of_rec(sigma_ref));
  return make_holder(
      [](auto&& sigma_ref, auto&& x, auto&& mu_ref) {
        return fma(std::forward<decltype(sigma_ref)>(sigma_ref),
                   std::forward<decltype(x)>(x),
                   std::forward<decltype(mu_ref)>(mu_ref));
      },
      std::forward<decltype(sigma_ref)>(sigma_ref), std::forward<T>(x),
      std::forward<decltype(mu_ref)>(mu_ref));
}

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * If the offset is zero and multiplier is one, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam T A scalar type or type inheriting from `Eigen::DenseBase`
 * @tparam M A scalar type or type inheriting from `Eigen::DenseBase`
 * @tparam S A scalar type or type inheriting from `Eigen::DenseBase`
 * @param[in] x Unconstrained scalar input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @param[in,out] lp Reference to log probability to increment.
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S, typename Lp,
          require_convertible_t<return_type_t<T, M, S>, Lp>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T, M, S>* = nullptr,
          require_all_not_std_vector_t<M, S>* = nullptr>
inline auto offset_multiplier_constrain(T&& x, M&& mu, S&& sigma, Lp& lp) {
  if constexpr (is_matrix_v<T> && is_matrix_v<M>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "mu", mu);
  }
  if constexpr (is_matrix_v<T> && is_matrix_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "sigma", sigma);
  }
  if constexpr (is_matrix_v<M> && is_matrix_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "mu", mu, "sigma",
                        sigma);
  }
  auto&& mu_ref = to_ref(std::forward<M>(mu));
  auto&& sigma_ref = to_ref(std::forward<S>(sigma));
  check_finite("offset_multiplier_constrain", "offset", value_of_rec(mu_ref));
  check_positive_finite("offset_multiplier_constrain", "multiplier",
                        value_of_rec(sigma_ref));
  if (stan::math::size(sigma_ref) == 1) {
    lp += sum(multiply_log(static_cast<double>(math::size(x)), sigma_ref));
  } else {
    lp += sum(log(sigma_ref));
  }
  return make_holder(
      [](auto&& sigma_ref, auto&& x, auto&& mu_ref) {
        return fma(std::forward<decltype(sigma_ref)>(sigma_ref),
                   std::forward<decltype(x)>(x),
                   std::forward<decltype(mu_ref)>(mu_ref));
      },
      std::forward<decltype(sigma_ref)>(sigma_ref), std::forward<T>(x),
      std::forward<decltype(mu_ref)>(mu_ref));
}

/**
 * Overload for when x and mu or sigma are `std::vectors`
 */
template <typename T, typename M, typename S,
          require_any_std_vector_t<T, M, S>* = nullptr>
inline auto offset_multiplier_constrain(const T& x, M&& mu, S&& sigma) {
  if constexpr (is_std_vector_v<T> && is_std_vector_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "sigma", sigma);
  }
  if constexpr (is_std_vector_v<T> && is_std_vector_v<M>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "mu", mu);
  }
  if constexpr (is_std_vector_v<M> && is_std_vector_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "mu", mu, "sigma",
                        sigma);
  }
  auto iter = [](auto&& it, std::size_t i) -> decltype(auto) {
    if constexpr (is_std_vector_v<decltype(it)>) {
      return it[i];
    } else {
      return it;
    }
  };
  auto&& mu_ref = to_ref(std::forward<M>(mu));
  auto&& sigma_ref = to_ref(std::forward<S>(sigma));
  using inner_ret_t = decltype(offset_multiplier_constrain(
      iter(x, 0), iter(mu_ref, 0), iter(sigma_ref, 0)));
  std::vector<plain_type_t<inner_ret_t>> ret;
  // In the language, if mu or sigma is a vector, x must also be a vector
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret.emplace_back(
        offset_multiplier_constrain(x[i], iter(mu_ref, i), iter(sigma_ref, i)));
  }
  return ret;
}

/**
 * Overload for when x and mu or sigma are `std::vectors`
 */
template <typename T, typename M, typename S, typename Lp,
          require_convertible_t<return_type_t<T, M, S>, Lp>* = nullptr,
          require_any_std_vector_t<T, M, S>* = nullptr>
inline auto offset_multiplier_constrain(const T& x, M&& mu, S&& sigma, Lp& lp) {
  if constexpr (is_std_vector_v<T> && is_std_vector_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "sigma", sigma);
  }
  if constexpr (is_std_vector_v<T> && is_std_vector_v<M>) {
    check_matching_dims("offset_multiplier_constrain", "x", x, "mu", mu);
  }
  if constexpr (is_std_vector_v<M> && is_std_vector_v<S>) {
    check_matching_dims("offset_multiplier_constrain", "mu", mu, "sigma",
                        sigma);
  }
  auto iter = [](auto&& it, std::size_t i) -> decltype(auto) {
    if constexpr (is_std_vector_v<decltype(it)>) {
      return it[i];
    } else {
      return it;
    }
  };
  auto&& mu_ref = to_ref(std::forward<M>(mu));
  auto&& sigma_ref = to_ref(std::forward<S>(sigma));
  using inner_ret_t = decltype(offset_multiplier_constrain(
      iter(x, 0), iter(mu_ref, 0), iter(sigma_ref, 0)));
  std::vector<plain_type_t<inner_ret_t>> ret;
  // In the language, if mu or sigma is a vector, x must also be a vector
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret.emplace_back(offset_multiplier_constrain(x[i], iter(mu_ref, i),
                                                 iter(sigma_ref, i), lp));
  }
  return ret;
}

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier. If the `Jacobian` parameter is `true`,
 * the log density accumulator is incremented with the log absolute Jacobian
 * determinant of the transform.  All of the transforms are specified with their
 * Jacobians in the *Stan Reference Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam M A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam S A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam Lp Scalar, the scalar types of T, M, and S should be convertable to
 * this
 * @param[in] x Unconstrained scalar input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @param[in, out] lp log density accumulator
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <bool Jacobian, typename T, typename M, typename S, typename Lp,
          require_convertible_t<return_type_t<T, M, S>, Lp>* = nullptr>
inline auto offset_multiplier_constrain(const T& x, const M& mu, const S& sigma,
                                        Lp& lp) {
  if constexpr (Jacobian) {
    return offset_multiplier_constrain(x, mu, sigma, lp);
  } else {
    return offset_multiplier_constrain(x, mu, sigma);
  }
}

}  // namespace math
}  // namespace stan

#endif
