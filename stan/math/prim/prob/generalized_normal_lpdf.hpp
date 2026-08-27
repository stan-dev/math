#ifndef STAN_MATH_PRIM_PROB_GENERALIZED_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GENERALIZED_NORMAL_LPDF_HPP

#include <cmath>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the generalized normal density for the specified scalar(s) given
 * the specified location, scale and shape parameters. y, mu, alpha, or beta can
 * each be either a scalar or a vector. Any vector inputs must be the same
 * length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/scale/shape tuple.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @tparam T_shape type of shape parameter
 * @param y (Sequence of) scalar(s)
 * @param mu (Sequence of) location parameter(s)
 * @param alpha (Sequence of) scale parameter(s)
 * @param beta (Sequence of) shape parameter(s)
 * @return The log of the product of the densities
 * @throw std::domain_error if alpha or beta is not positive or not finite
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_shape> generalized_normal_lpdf(
    T_y&& y, T_loc&& mu, T_scale&& alpha, T_shape&& beta) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_scale>;
  using T_beta_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "generalized_normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", alpha, "Shape parameter", beta);

  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_alpha_ref alpha_ref = std::forward<T_scale>(alpha);
  T_beta_ref beta_ref = std::forward<T_shape>(beta);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", alpha_val);

  // With β = +∞ this could be defined to be uniform, but we don't support that.
  check_positive_finite(function, "Shape parameter", beta_val);

  if (size_zero(y, mu, alpha, beta)) {
    return 0;
  }
  if constexpr (!include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    return 0;
  }

  const auto& one_plus_inv_beta
      = to_ref_if<!is_constant<T_shape>::value>(inv(beta_val) + 1);
  const auto& residual
      = to_ref_if<!is_constant_all<T_y, T_loc>::value>(y_val - mu_val);
  const auto& inverse_scale
      = to_ref_if<!is_constant_all<T_y, T_loc>::value>(inv(alpha_val));
  const auto& scaled_residual = to_ref(residual / alpha_val);
  const auto& abs_scaled_residual = to_ref(abs(scaled_residual));
  const auto& scaled_residual_is_zero
      = eval(value_of_rec(abs_scaled_residual) == 0);
  auto abs_scaled_residual_or_one = eval(abs_scaled_residual);
  if constexpr (is_eigen_v<std::decay_t<decltype(scaled_residual_is_zero)>>) {
    abs_scaled_residual_or_one
        = scaled_residual_is_zero.select(1.0, abs_scaled_residual_or_one);
  } else {
    if (scaled_residual_is_zero) {
      abs_scaled_residual_or_one += 1;
    }
  }
  auto abs_scaled_residual_pow_beta
      = eval(pow(abs_scaled_residual_or_one, beta_val));
  // At beta == 2 the kernel is exactly square(scaled_residual). Retaining that
  // expression at the tie preserves the normal Hessian; differentiating the
  // generic abs/pow expression would instead inherit abs'(0) == 0.
  auto beta_is_two = eval(value_of_rec(beta_val) == 2.0);
  auto apply_zero_residual_branch = [&](auto& generic_value,
                                        const auto& beta_two_value) {
    if constexpr (is_eigen_v<std::decay_t<decltype(scaled_residual_is_zero)>>) {
      if constexpr (is_eigen_v<std::decay_t<decltype(beta_is_two)>>) {
        generic_value = scaled_residual_is_zero.select(
            beta_is_two.select(beta_two_value, 0.0), generic_value);
      } else if (beta_is_two) {
        generic_value
            = scaled_residual_is_zero.select(beta_two_value, generic_value);
      } else {
        generic_value = scaled_residual_is_zero.select(0.0, generic_value);
      }
    } else if (scaled_residual_is_zero) {
      if constexpr (is_eigen_v<std::decay_t<decltype(beta_is_two)>>) {
        generic_value = beta_is_two.select(beta_two_value, 0.0);
      } else {
        generic_value = beta_is_two ? beta_two_value : 0;
      }
    }
  };
  auto beta_two_kernel = eval(square(scaled_residual) + 0.0 * beta_val);
  apply_zero_residual_branch(abs_scaled_residual_pow_beta, beta_two_kernel);
  const size_t num_terms = max_size(y, mu, alpha, beta);

  T_partials_return logp = -sum(abs_scaled_residual_pow_beta);

  if constexpr (include_summand<propto>::value) {
    logp -= LOG_TWO * num_terms;
  }
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= sum(log(alpha_val)) * (num_terms / math::size(alpha));
  }
  if constexpr (include_summand<propto, T_shape>::value) {
    logp -= sum(lgamma(one_plus_inv_beta)) * (num_terms / math::size(beta));
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, alpha_ref, beta_ref);

  if constexpr (!is_constant_all<T_y, T_loc>::value) {
    // At y == mu, the derivative is 0 for beta > 1 and undefined otherwise.
    // Use abs_scaled_residual_or_one (0 replaced with 1) to avoid NaN from
    // pow(0, beta-1), and preserve the exact normal curvature when beta == 2.
    auto kernel_derivative_wrt_residual
        = eval(sign(residual) * beta_val
               * pow(abs_scaled_residual_or_one, beta_val - 1) * inverse_scale);
    const auto& beta_two_kernel_derivative
        = eval(2.0 * scaled_residual / alpha_val + 0.0 * beta_val);
    apply_zero_residual_branch(kernel_derivative_wrt_residual,
                               beta_two_kernel_derivative);
    if constexpr (!is_constant<T_y>::value) {
      partials<0>(ops_partials) = -kernel_derivative_wrt_residual;
    }
    if constexpr (!is_constant<T_loc>::value) {
      partials<1>(ops_partials) = std::move(kernel_derivative_wrt_residual);
    }
  }
  if constexpr (!is_constant<T_scale>::value) {
    partials<2>(ops_partials)
        = (beta_val * abs_scaled_residual_pow_beta - 1) * inverse_scale;
  }
  if constexpr (!is_constant<T_shape>::value) {
    // multiply_log(0, 0) = 0 by convention, but fvar autodiff of
    // multiply_log at (0, 0) produces NaN. abs_scaled_residual_or_one avoids
    // this.
    partials<3>(ops_partials)
        = digamma(one_plus_inv_beta) * inv_square(beta_val)
          - multiply_log(abs_scaled_residual_pow_beta,
                         abs_scaled_residual_or_one);
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> generalized_normal_lpdf(
    T_y&& y, T_loc&& mu, T_scale&& alpha, T_shape&& beta) {
  return generalized_normal_lpdf<false>(
      std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(alpha), std::forward<T_shape>(beta));
}

}  // namespace math
}  // namespace stan
#endif
