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
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
// Exposing to let me use in tests
// The current tests fail for 1e8 and pass for 1e9, so setting to 1e10
constexpr double neg_binomial_alpha_cutoff = 1e10;
}  // namespace internal

// NegBinomial(n|alpha, beta)  [alpha > 0;  beta > 0;  n >= 0]
template <bool propto, typename T_n, typename T_shape, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_shape, T_inv_scale>* = nullptr>
return_type_t<T_shape, T_inv_scale> neg_binomial_lpmf(const T_n& n,
                                                      const T_shape& alpha,
                                                      const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;
  using std::log;
  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  static const char* function = "neg_binomial_lpmf";
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "Failures variable", n_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  if (size_zero(n, alpha, beta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_shape, T_inv_scale>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  operands_and_partials<T_alpha_ref, T_beta_ref> ops_partials(alpha_ref,
                                                              beta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_alpha = stan::math::size(alpha);
  size_t size_beta = stan::math::size(beta);
  size_t size_alpha_beta = max_size(alpha, beta);
  size_t max_size_seq_view = max_size(n, alpha, beta);

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_alpha(size_alpha);
  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_alpha; ++i) {
      digamma_alpha[i] = digamma(alpha_vec.val(i));
    }
  }

  VectorBuilder<true, T_partials_return, T_inv_scale> log1p_inv_beta(size_beta);
  VectorBuilder<true, T_partials_return, T_inv_scale> log1p_beta(size_beta);
  for (size_t i = 0; i < size_beta; ++i) {
    const T_partials_return beta_dbl = beta_vec.val(i);
    log1p_inv_beta[i] = log1p(inv(beta_dbl));
    log1p_beta[i] = log1p(beta_dbl);
  }

  VectorBuilder<!is_constant_all<T_inv_scale>::value, T_partials_return,
                T_shape, T_inv_scale>
      lambda_m_alpha_over_1p_beta(size_alpha_beta);
  if (!is_constant_all<T_inv_scale>::value) {
    for (size_t i = 0; i < size_alpha_beta; ++i) {
      const T_partials_return alpha_dbl = alpha_vec.val(i);
      const T_partials_return beta_dbl = beta_vec.val(i);
      lambda_m_alpha_over_1p_beta[i]
          = alpha_dbl / beta_dbl - alpha_dbl / (1 + beta_dbl);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const T_partials_return alpha_dbl = alpha_vec.val(i);
    const T_partials_return beta_dbl = beta_vec.val(i);

    if (include_summand<propto, T_shape>::value) {
      if (n_vec[i] != 0) {
        logp += binomial_coefficient_log(n_vec[i] + alpha_dbl - 1.0,
                                         alpha_dbl - 1.0);
      }
    }
    logp -= alpha_dbl * log1p_inv_beta[i] + n_vec[i] * log1p_beta[i];

    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge1_.partials_[i] += digamma(alpha_dbl + n_vec[i])
                                          - digamma_alpha[i]
                                          - log1p_inv_beta[i];
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[i]
          += lambda_m_alpha_over_1p_beta[i] - n_vec[i] / (beta_dbl + 1.0);
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
