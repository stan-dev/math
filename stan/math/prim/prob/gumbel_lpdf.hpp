#ifndef STAN_MATH_PRIM_PROB_GUMBEL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GUMBEL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
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
template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> gumbel_lpdf(const T_y& y, const T_loc& mu,
                                               const T_scale& beta) {
  static const char* function = "gumbel_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;

  using std::exp;
  using std::log;

  if (size_zero(y, mu, beta)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive(function, "Scale parameter", beta);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", beta);

  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> beta_vec(beta);
  size_t N = max_size(y, mu, beta);

  VectorBuilder<true, T_partials_return, T_scale> inv_beta(size(beta));
  VectorBuilder<include_summand<propto, T_scale>::value, T_partials_return,
                T_scale>
      log_beta(size(beta));
  for (size_t i = 0; i < size(beta); i++) {
    inv_beta[i] = 1.0 / value_of(beta_vec[i]);
    if (include_summand<propto, T_scale>::value) {
      log_beta[i] = log(value_of(beta_vec[i]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);

    const T_partials_return y_minus_mu_over_beta
        = (y_dbl - mu_dbl) * inv_beta[n];

    if (include_summand<propto, T_scale>::value) {
      logp -= log_beta[n];
    }
    logp += -y_minus_mu_over_beta - exp(-y_minus_mu_over_beta);

    T_partials_return scaled_diff = inv_beta[n] * exp(-y_minus_mu_over_beta);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= inv_beta[n] - scaled_diff;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += inv_beta[n] - scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] += -inv_beta[n]
                                          + y_minus_mu_over_beta * inv_beta[n]
                                          - scaled_diff * y_minus_mu_over_beta;
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
