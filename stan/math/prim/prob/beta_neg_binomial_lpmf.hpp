#ifndef STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_rising_factorial_ratio.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the Beta Negative Binomial distribution with given
 * number of successes, prior success, and prior failure parameters.
 * Given containers of matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_n type of failure parameter
 * @tparam T_r type of number of successes parameter
 * @tparam T_alpha type of prior success parameter
 * @tparam T_beta type of prior failure parameter
 *
 * @param n failure parameter
 * @param r Number of successes parameter
 * @param alpha prior success parameter
 * @param beta prior failure parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if r, alpha, or beta fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_r, typename T_alpha,
          typename T_beta,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_r, T_alpha, T_beta>* = nullptr>
inline return_type_t<T_r, T_alpha, T_beta> beta_neg_binomial_lpmf(
    const T_n& n, const T_r& r, const T_alpha& alpha, const T_beta& beta) {
  using T_partials_return = partials_return_t<T_n, T_r, T_alpha, T_beta>;
  using T_n_ref = ref_type_t<T_n>;
  using T_r_ref = ref_type_t<T_r>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_beta_ref = ref_type_t<T_beta>;
  static constexpr const char* function = "beta_neg_binomial_lpmf";
  check_consistent_sizes(
      function, "Failures variable", n, "Number of successes parameter", r,
      "Prior success parameter", alpha, "Prior failure parameter", beta);
  if (size_zero(n, r, alpha, beta)) {
    return 0.0;
  }

  T_n_ref n_ref = n;
  T_r_ref r_ref = r;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "Failures variable", n_ref);
  check_positive_finite(function, "Number of successes parameter", r_ref);
  check_positive_finite(function, "Prior success parameter", alpha_ref);
  check_positive_finite(function, "Prior failure parameter", beta_ref);

  if constexpr (!include_summand<propto, T_r, T_alpha, T_beta>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(r_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_r_ref> r_vec(r_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  const size_t max_size_seq_view = max_size(n, r, alpha, beta);
  // With D(x, k) = lgamma(x + k) - lgamma(x) and A = alpha + beta, the log
  // pmf is
  //
  //   -lgamma(n + 1) + D(beta, n) + [D(r, n) - D(r + A, n)]
  //   + [D(alpha, r) - D(A, r)].
  //
  // Formed from lbeta and lgamma directly, the terms are of the size of the
  // shapes and cancel when the shapes are large. Here
  // D(beta, n) = lgamma(n) - lbeta(n, beta) for n > 0, and each bracket is
  // internal::log_rising_factorial_ratio, which does not cancel. The
  // partials are sums of digamma differences psi(x + k) - psi(x), formed by
  // digamma_diff.
  T_partials_return logp(0.0);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    const T_partials_return r_dbl = r_vec.val(i);
    const T_partials_return alpha_dbl = alpha_vec.val(i);
    const T_partials_return beta_dbl = beta_vec.val(i);
    const T_partials_return n_dbl = n_vec.val(i);
    const T_partials_return alpha_plus_beta = alpha_dbl + beta_dbl;
    if (n_dbl > 0) {
      if constexpr (include_summand<propto>::value) {
        logp -= log(n_dbl);  // -lgamma(n + 1) + lgamma(n)
      } else {
        logp += lgamma(n_dbl);
      }
      logp -= lbeta(n_dbl, beta_dbl);
    }
    logp += internal::log_rising_factorial_ratio(r_dbl, alpha_plus_beta, n_dbl)
            + internal::log_rising_factorial_ratio(alpha_dbl, beta_dbl, r_dbl);

    if constexpr (is_any_autodiff_v<T_r, T_alpha, T_beta>) {
      const T_partials_return dpsi_r_plus_ab_n
          = digamma_diff(r_dbl + alpha_plus_beta, n_dbl);
      if constexpr (is_autodiff_v<T_r>) {
        partials<0>(ops_partials)[i]
            += digamma_diff(r_dbl, n_dbl) - dpsi_r_plus_ab_n
               - digamma_diff(r_dbl + alpha_dbl, beta_dbl);
      }
      if constexpr (is_any_autodiff_v<T_alpha, T_beta>) {
        const T_partials_return dpsi_ab_r
            = digamma_diff(alpha_plus_beta, r_dbl);
        if constexpr (is_autodiff_v<T_alpha>) {
          partials<1>(ops_partials)[i]
              += digamma_diff(alpha_dbl, r_dbl) - dpsi_ab_r - dpsi_r_plus_ab_n;
        }
        if constexpr (is_autodiff_v<T_beta>) {
          partials<2>(ops_partials)[i]
              += digamma_diff(beta_dbl, n_dbl) - dpsi_r_plus_ab_n - dpsi_ab_r;
        }
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_r, typename T_alpha, typename T_beta>
inline return_type_t<T_r, T_alpha, T_beta> beta_neg_binomial_lpmf(
    const T_n& n, const T_r& r, const T_alpha& alpha, const T_beta& beta) {
  return beta_neg_binomial_lpmf<false>(n, r, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
