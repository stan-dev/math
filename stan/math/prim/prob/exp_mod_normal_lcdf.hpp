#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>

namespace stan {
namespace math {
namespace internal {

/** log F = log_diff_exp(a, b) with a = log Phi(z), b = v^2/2 - lambda (y - mu)
 * + log Phi(z - v), v = lambda sigma; log (1 - F) = log_sum_exp(log Phi(-z),
 * b).
 */
template <bool upper, typename T_y, typename T_loc, typename T_scale,
          typename T_inv_scale>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lcdf_impl(
    const char* function, T_y&& y, T_loc&& mu, T_scale&& sigma,
    T_inv_scale&& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  constexpr bool any_autodiff
      = is_any_autodiff_v<T_y, T_loc, T_scale, T_inv_scale>;
  constexpr double sign = upper ? -1.0 : 1.0;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale parameter",
                         lambda);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  T_lambda_ref lambda_ref = std::forward<T_inv_scale>(lambda);
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));
  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);
  check_positive_finite(function, "Inv_scale parameter", lambda_val);

  if (size_zero(y_ref, mu_ref, sigma_ref, lambda_ref)) {
    return 0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);
  if (any(y_val == NEGATIVE_INFTY)) {
    return ops_partials.build(upper ? 0.0 : NEGATIVE_INFTY);
  }
  if (any(y_val == INFTY)) {
    return ops_partials.build(upper ? NEGATIVE_INFTY : 0.0);
  }

  const auto& diff = to_ref(y_val - mu_val);
  const auto& z = to_ref(diff / sigma_val);
  const auto& v = to_ref(lambda_val * sigma_val);
  const auto [log_a, slope_a]
      = internal::std_normal_lcdf_value_grad<any_autodiff>(sign * z);
  const auto [log_phi_b, slope_b]
      = internal::std_normal_lcdf_value_grad<any_autodiff>(z - v);
  const auto& log_b = to_ref_if<any_autodiff>(0.5 * square(v)
                                              - lambda_val * diff + log_phi_b);
  const auto& lp = to_ref_if<any_autodiff>([&]() {
    if constexpr (upper) {
      return log_sum_exp(log_a, log_b);
    } else {
      return log_diff_exp(log_a, log_b);
    }
  }());
  const T_partials_return cdf_log = sum(lp);

  if constexpr (any_autodiff) {
    // Weights of the two terms in the total; the second is signed.
    const auto& w_b = to_ref(-sign * exp(log_b - lp));
    const auto& s_b = to_ref(w_b * slope_b);
    const auto& s
        = to_ref_if<(is_autodiff_v<T_y> + is_autodiff_v<T_loc>
                     + is_autodiff_v<T_scale>)
                    >= 2>((sign * exp(log_a - lp) * slope_a + s_b) / sigma_val);
    const auto& q
        = to_ref_if<is_all_autodiff_v<T_scale, T_inv_scale>>(w_b * v - s_b);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = s - w_b * lambda_val;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = w_b * lambda_val - s;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      partials<2>(ops_partials) = select(s == 0, 0.0, -s * z) + lambda_val * q;
    }
    if constexpr (is_autodiff_v<T_inv_scale>) {
      partials<3>(ops_partials) = sigma_val * q - w_b * diff;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace internal

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lcdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_inv_scale&& lambda) {
  return internal::exp_mod_normal_lcdf_impl<false>(
      "exp_mod_normal_lcdf", std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma), std::forward<T_inv_scale>(lambda));
}

}  // namespace math
}  // namespace stan
#endif
