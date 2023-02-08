#ifndef STAN_MATH_PRIM_PROB_GUMBEL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GUMBEL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the Gumbel log probability density for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of densities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param mu location parameter
 * @param beta scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is nan, mu is infinite, or beta is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> gumbel_lpdf(const T_y& y, const T_loc& mu,
                                               const T_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "gumbel_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", beta);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", beta_val);

  if (size_zero(y, mu, beta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  operands_and_partials<T_y_ref, T_mu_ref, T_beta_ref> ops_partials(
      y_ref, mu_ref, beta_ref);

  const auto& inv_beta
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(inv(beta_val));
  const auto& y_minus_mu_over_beta = to_ref((y_val - mu_val) * inv_beta);
  const auto& exp_y_m_mu_over_beta
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
          exp(-y_minus_mu_over_beta));

  size_t N = max_size(y, mu, beta);
  T_partials_return logp = -sum(y_minus_mu_over_beta + exp_y_m_mu_over_beta);
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(beta_val)) * N / math::size(beta);
  }

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    const auto& scaled_diff
        = to_ref_if<!is_constant_all<T_loc>::value
                        + !is_constant_all<T_scale>::value
                        + !is_constant_all<T_y>::value
                    >= 2>(inv_beta * exp_y_m_mu_over_beta - inv_beta);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = scaled_diff;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_ = -scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_
          = -y_minus_mu_over_beta * scaled_diff - inv_beta;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> gumbel_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& beta) {
  return gumbel_lpdf<false>(y, mu, beta);
}

}  // namespace math
}  // namespace stan
#endif
