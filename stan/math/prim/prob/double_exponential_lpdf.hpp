#ifndef STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the double exponential log probability density function. Given
 * containers of matching sizes, returns the log sum of densities.
 *
 * @tparam T_y type of real parameter.
 * @tparam T_loc type of location parameter.
 * @tparam T_scale type of scale parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is nan, mu is infinite, or sigma is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> double_exponential_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "double_exponential_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Shape parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_finite(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  const auto& inv_sigma = to_ref(inv(sigma_val));
  const auto& y_m_mu
      = to_ref_if<!is_constant_all<T_y, T_loc>::value>(y_val - mu_val);
  const auto& abs_diff_y_mu = fabs(y_m_mu);
  const auto& scaled_diff
      = to_ref_if<!is_constant_all<T_scale>::value>(abs_diff_y_mu * inv_sigma);

  size_t N = max_size(y, mu, sigma);
  if (include_summand<propto>::value) {
    logp -= N * LOG_TWO;
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }
  logp -= sum(scaled_diff);

  if (!is_constant_all<T_y, T_loc>::value) {
    const auto& diff_sign = sign(y_m_mu);
    const auto& rep_deriv
        = to_ref_if<(!is_constant_all<T_y>::value
                     && !is_constant_all<T_loc>::value)>(diff_sign * inv_sigma);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = -rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) = rep_deriv;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    partials<2>(ops_partials) = inv_sigma * (scaled_diff - 1);
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> double_exponential_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  return double_exponential_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
