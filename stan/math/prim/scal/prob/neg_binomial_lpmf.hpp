#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

// NegBinomial(n|alpha, beta)  [alpha > 0;  beta > 0;  n >= 0]
template <bool propto, typename T_n, typename T_shape, typename T_inv_scale>
inline auto neg_binomial_lpmf(T_n&& n, T_shape&& alpha,
                              T_inv_scale&& beta) {
  using T_partials = partials_return_t<T_n, T_shape, T_inv_scale>;
  T_partials logp(0.0);

  static const char* function = "neg_binomial_lpmf";
  check_nonnegative(function, "Failures variable", n);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  using std::log;

  const scalar_seq_view<T_n> n_vec(n);
  const scalar_seq_view<T_shape> alpha_vec(alpha);
  const scalar_seq_view<T_inv_scale> beta_vec(beta);
  const size_t size = max_size(n, alpha, beta);

  operands_and_partials<T_shape, T_inv_scale> ops_partials(alpha, beta);
  if (!include_summand<propto, T_shape, T_inv_scale>::value) {
    return ops_partials.build(logp);
  } else if (size_zero(n, alpha, beta)) {
    return ops_partials.build(logp);
  }

  for (size_t i = 0; i < size; i++) {
    auto lambda = value_of(alpha_vec[i]) / value_of(beta_vec[i]);
    if (alpha_vec[i] > 1e10) {  // reduces numerically to Poisson
      if (include_summand<propto>::value) {
        logp -= lgamma(n_vec[i] + 1.0);
      }
      if (include_summand<propto, T_shape, T_inv_scale>::value) {
        logp += multiply_log(n_vec[i], lambda) - lambda;
      }

      if (!is_constant_all<T_shape>::value) {
        ops_partials.edge1_.partials_[i]
            += n_vec[i] / value_of(alpha_vec[i]) - 1.0 / value_of(beta_vec[i]);
      }
      if (!is_constant_all<T_inv_scale>::value) {
        ops_partials.edge2_.partials_[i]
            += (lambda - n_vec[i]) / value_of(beta_vec[i]);
      }
    } else {  // standard density definition
      if (include_summand<propto, T_shape>::value) {
        if (n_vec[i] != 0) {
          logp += binomial_coefficient_log(
              n_vec[i] + value_of(alpha_vec[i]) - 1.0, n_vec[i]);
        }
      }
      auto log1p_beta = log1p(value_of(beta_vec[i]));

      if (include_summand<propto, T_shape, T_inv_scale>::value) {
        auto alpha_times_log_beta_over_1p_beta
            = value_of(alpha_vec[i])
              * log(value_of(beta_vec[i]) / (1.0 + value_of(beta_vec[i])));
        logp += alpha_times_log_beta_over_1p_beta - n_vec[i] * log1p_beta;
      }

      if (!is_constant_all<T_shape>::value) {
        auto digamma_alpha = digamma(value_of(alpha_vec[i]));
        auto log_beta_m_log1p_beta = log(value_of(beta_vec[i])) - log1p_beta;
        ops_partials.edge1_.partials_[i]
            += digamma(value_of(alpha_vec[i]) + n_vec[i]) - digamma_alpha
               + log_beta_m_log1p_beta;
      }
      if (!is_constant_all<T_inv_scale>::value) {
        auto lambda_m_alpha_over_1p_beta
            = lambda
              - (value_of(alpha_vec[i]) / (1.0 + value_of(beta_vec[i])));
        ops_partials.edge2_.partials_[i]
            += lambda_m_alpha_over_1p_beta
               - n_vec[i] / (value_of(beta_vec[i]) + 1.0);
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_shape, typename T_inv_scale>
inline auto neg_binomial_lpmf(T_n&& n, T_shape&& alpha,
                              T_inv_scale&& beta) {
  return neg_binomial_lpmf<false>(std::forward<T_n>(n), std::forward<T_shape>(alpha), std::forward<T_inv_scale>(beta));
}

}  // namespace math
}  // namespace stan
#endif
