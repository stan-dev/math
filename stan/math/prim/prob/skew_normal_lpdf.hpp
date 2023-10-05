#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  static const char* function = "skew_normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape paramter", alpha);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_alpha_ref alpha_ref = alpha;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_finite(function, "Shape parameter", alpha_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma, alpha)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    return 0.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, alpha_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(inv(sigma_val));
  const auto& y_minus_mu_over_sigma
      = to_ref_if<include_summand<propto, T_y, T_loc, T_scale, T_shape>::value>(
          (y_val - mu_val) * inv_sigma);
  const auto& log_erfc_alpha_z
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale, T_shape>::value>(
          log(erfc(-alpha_val * y_minus_mu_over_sigma * INV_SQRT_TWO)));

  size_t N = max_size(y, mu, sigma, alpha);
  T_partials_return logp = sum(log_erfc_alpha_z);
  if (include_summand<propto>::value) {
    logp -= HALF_LOG_TWO_PI * N;
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }
  if (include_summand<propto, T_y, T_loc, T_scale>::value) {
    logp -= sum(square(y_minus_mu_over_sigma)) * 0.5 * N
            / max_size(y, mu, sigma);
  }

  if (!is_constant_all<T_y, T_loc, T_scale, T_shape>::value) {
    const auto& sq = square(alpha_val * y_minus_mu_over_sigma * INV_SQRT_TWO);
    const auto& ex = exp(-sq - log_erfc_alpha_z);
    auto deriv_logerf = to_ref_if<!is_constant_all<T_y, T_loc>::value
                                      + !is_constant_all<T_scale>::value
                                      + !is_constant_all<T_shape>::value
                                  >= 2>(SQRT_TWO_OVER_SQRT_PI * ex);
    if (!is_constant_all<T_y, T_loc>::value) {
      auto deriv_y_loc = to_ref_if<(!is_constant_all<T_y>::value
                                    && !is_constant_all<T_loc>::value)>(
          (y_minus_mu_over_sigma - deriv_logerf * alpha_val) * inv_sigma);
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials) = -deriv_y_loc;
      }
      if (!is_constant_all<T_loc>::value) {
        partials<1>(ops_partials) = std::move(deriv_y_loc);
      }
    }
    if (!is_constant_all<T_scale>::value) {
      edge<2>(ops_partials).partials_
          = ((y_minus_mu_over_sigma - deriv_logerf * alpha_val)
                 * y_minus_mu_over_sigma
             - 1)
            * inv_sigma;
    }
    if (!is_constant_all<T_shape>::value) {
      partials<3>(ops_partials) = deriv_logerf * y_minus_mu_over_sigma;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  return skew_normal_lpdf<false>(y, mu, sigma, alpha);
}

}  // namespace math
}  // namespace stan
#endif
