#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
// Exposing to let me use in tests
// The current tests fail for 1e8 and pass for 1e9, so setting to 1e10
constexpr double neg_binomial_alpha_cutoff = 1e10;
}  // namespace internal

// NegBinomial(n|alpha, beta)  [alpha > 0;  beta > 0;  n >= 0]
template <bool propto, typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> neg_binomial_lpmf(const T_n& n,
                                                      const T_shape& alpha,
                                                      const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;

  static const char* function = "neg_binomial_lpmf";

  if (size_zero(n, alpha, beta)) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  check_nonnegative(function, "Failures variable", n);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  if (!include_summand<propto, T_shape, T_inv_scale>::value) {
    return 0.0;
  }

  using std::log;

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t size_alpha = size(alpha);
  size_t size_beta = size(beta);
  size_t size_alpha_beta = max_size(alpha, beta);
  size_t max_size_seq_view = max_size(n, alpha, beta);

  operands_and_partials<T_shape, T_inv_scale> ops_partials(alpha, beta);

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_alpha(size_alpha);
  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_alpha; ++i) {
      digamma_alpha[i] = digamma(value_of(alpha_vec[i]));
    }
  }

  VectorBuilder<true, T_partials_return, T_inv_scale> log_beta(size_beta);
  VectorBuilder<true, T_partials_return, T_inv_scale> log1p_beta(size_beta);
  VectorBuilder<true, T_partials_return, T_inv_scale> log_beta_m_log1p_beta(
      size_beta);
  for (size_t i = 0; i < size_beta; ++i) {
    const T_partials_return beta_dbl = value_of(beta_vec[i]);
    log_beta[i] = log(beta_dbl);
    log1p_beta[i] = log1p(beta_dbl);
    log_beta_m_log1p_beta[i] = log_beta[i] - log1p_beta[i];
  }

  VectorBuilder<true, T_partials_return, T_shape, T_inv_scale> lambda(
      size_alpha_beta);
  VectorBuilder<true, T_partials_return, T_shape, T_inv_scale>
      alpha_log_beta_over_1p_beta(size_alpha_beta);
  VectorBuilder<!is_constant_all<T_inv_scale>::value, T_partials_return,
                T_shape, T_inv_scale>
      lambda_m_alpha_over_1p_beta(size_alpha_beta);
  for (size_t i = 0; i < size_alpha_beta; ++i) {
    const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
    const T_partials_return beta_dbl = value_of(beta_vec[i]);
    lambda[i] = alpha_dbl / beta_dbl;
    alpha_log_beta_over_1p_beta[i] = alpha_dbl * log_beta_m_log1p_beta[i];
    if (!is_constant_all<T_inv_scale>::value) {
      lambda_m_alpha_over_1p_beta[i] = lambda[i] - alpha_dbl / (1 + beta_dbl);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    if (alpha_vec[i] > internal::neg_binomial_alpha_cutoff) {
      // reduces numerically to Poisson
      if (include_summand<propto>::value) {
        logp -= lgamma(n_vec[i] + 1.0);
      }
      logp += multiply_log(n_vec[i], lambda[i]) - lambda[i];

      if (!is_constant_all<T_shape>::value) {
        ops_partials.edge1_.partials_[i]
            += n_vec[i] / value_of(alpha_vec[i]) - 1.0 / value_of(beta_vec[i]);
      }
      if (!is_constant_all<T_inv_scale>::value) {
        ops_partials.edge2_.partials_[i]
            += (lambda[i] - n_vec[i]) / value_of(beta_vec[i]);
      }
    } else {  // standard density definition
      if (include_summand<propto, T_shape>::value) {
        if (n_vec[i] != 0) {
          logp += binomial_coefficient_log(
              n_vec[i] + value_of(alpha_vec[i]) - 1.0, n_vec[i]);
        }
      }
      logp += alpha_log_beta_over_1p_beta[i] - n_vec[i] * log1p_beta[i];

      if (!is_constant_all<T_shape>::value) {
        ops_partials.edge1_.partials_[i]
            += digamma(value_of(alpha_vec[i]) + n_vec[i]) - digamma_alpha[i]
               + log_beta_m_log1p_beta[i];
      }
      if (!is_constant_all<T_inv_scale>::value) {
        ops_partials.edge2_.partials_[i]
            += lambda_m_alpha_over_1p_beta[i]
               - n_vec[i] / (value_of(beta_vec[i]) + 1.0);
      }
    }
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_shape, typename T_inv_scale>
inline return_type_t<T_shape, T_inv_scale> neg_binomial_lpmf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  return neg_binomial_lpmf<false>(n, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
