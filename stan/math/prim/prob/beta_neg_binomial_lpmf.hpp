#ifndef STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
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
  T_partials_return logp(0.0);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    if constexpr (include_summand<propto>::value) {
      logp -= lgamma(n_vec[i] + 1);
    }
    T_partials_return lbeta_denominator = lbeta(r_vec.val(i), alpha_vec.val(i));
    T_partials_return lgamma_numerator = lgamma(n_vec[i] + beta_vec.val(i));
    T_partials_return lgamma_denominator = lgamma(beta_vec.val(i));
    T_partials_return lbeta_numerator
        = lbeta(n_vec[i] + r_vec.val(i), alpha_vec.val(i) + beta_vec.val(i));
    logp += lbeta_numerator + lgamma_numerator - lbeta_denominator
            - lgamma_denominator;
    if (!is_constant_all<T_r, T_alpha, T_beta>::value) {
      T_partials_return digamma_n_r_alpha_beta = digamma(
          n_vec[i] + r_vec.val(i) + alpha_vec.val(i) + beta_vec.val(i));

      if constexpr (!is_constant<T_r>::value || !is_constant<T_alpha>::value) {
        T_partials_return digamma_r_alpha
            = digamma(r_vec.val(i) + alpha_vec.val(i));
        if constexpr (!is_constant_all<T_r>::value) {
          partials<0>(ops_partials)[i]
              += digamma(n_vec[i] + r_vec.val(i)) - digamma_n_r_alpha_beta
                 - (digamma(r_vec.val(i)) - digamma_r_alpha);
        }
        if constexpr (!is_constant_all<T_alpha>::value) {
          partials<1>(ops_partials)[i]
              += -digamma_n_r_alpha_beta
                 - (digamma(alpha_vec.val(i)) - digamma_r_alpha);
        }
      }
      if constexpr (!is_constant<T_beta>::value
                    || !is_constant<T_alpha>::value) {
        T_partials_return digamma_alpha_beta
            = digamma(alpha_vec.val(i) + beta_vec.val(i));
        if constexpr (!is_constant_all<T_beta>::value) {
          partials<2>(ops_partials)[i] += digamma_alpha_beta
                                          - digamma_n_r_alpha_beta
                                          + digamma(n_vec[i] + beta_vec.val(i))
                                          - digamma(beta_vec.val(i));
        }
        if constexpr (!is_constant_all<T_alpha>::value) {
          partials<1>(ops_partials)[i] += digamma_alpha_beta;
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
