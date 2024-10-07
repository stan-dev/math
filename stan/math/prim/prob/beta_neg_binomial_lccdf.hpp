#ifndef STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/hypergeometric_3F2.hpp>
#include <stan/math/prim/fun/grad_F32.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the Beta-Negative Binomial distribution with given
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
template <typename T_n, typename T_r, typename T_alpha, typename T_beta>
inline return_type_t<T_r, T_alpha, T_beta> beta_neg_binomial_lccdf(
    const T_n& n, const T_r& r, const T_alpha& alpha, const T_beta& beta) {
  using std::exp;
  using std::log;
  using T_partials_return = partials_return_t<T_n, T_r, T_alpha, T_beta>;
  using T_r_ref = ref_type_t<T_r>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_beta_ref = ref_type_t<T_beta>;
  static constexpr const char* function = "beta_neg_binomial_lccdf";
  check_consistent_sizes(
      function, "Failures variable", n, "Number of successes parameter", r,
      "Prior success parameter", alpha, "Prior failure parameter", beta);
  if (size_zero(n, r, alpha, beta)) {
    return 0;
  }

  T_r_ref r_ref = r;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Number of successes parameter", r_ref);
  check_positive_finite(function, "Prior success parameter", alpha_ref);
  check_positive_finite(function, "Prior failure parameter", beta_ref);

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(r_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_r_ref> r_vec(r_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_n = stan::math::size(n);
  size_t max_size_seq_view = max_size(n, r, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_n; i++) {
    if (n_vec.val(i) < 0) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit return for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(negative_infinity());
    }
    T_partials_return n_dbl = n_vec.val(i);
    T_partials_return r_dbl = r_vec.val(i);
    T_partials_return alpha_dbl = alpha_vec.val(i);
    T_partials_return beta_dbl = beta_vec.val(i);
    T_partials_return b_plus_n = beta_dbl + n_dbl;
    T_partials_return r_plus_n = r_dbl + n_dbl;
    T_partials_return a_plus_r = alpha_dbl + r_dbl;
    T_partials_return one = 1;
    T_partials_return precision = 1e-8;  // default -6, set -8 to pass all tests

    T_partials_return F
        = hypergeometric_3F2({one, b_plus_n + 1, r_plus_n + 1},
                             {n_dbl + 2, a_plus_r + b_plus_n + 1}, one);
    T_partials_return C = lgamma(r_plus_n + 1) + lbeta(a_plus_r, b_plus_n + 1)
                          - lgamma(r_dbl) - lbeta(alpha_dbl, beta_dbl)
                          - lgamma(n_dbl + 2);
    T_partials_return ccdf = exp(C) * F;
    T_partials_return P_i = log(ccdf);
    P += P_i;

    if (!is_constant_all<T_r, T_alpha, T_beta>::value) {
      T_partials_return digamma_n_r_alpha_beta
          = digamma(a_plus_r + b_plus_n + 1);
      T_partials_return dF[6];
      grad_F32(dF, one, b_plus_n + 1, r_plus_n + 1, n_dbl + 2,
               a_plus_r + b_plus_n + 1, one, precision, 1e5);

      if constexpr (!is_constant<T_r>::value || !is_constant<T_alpha>::value) {
        T_partials_return digamma_r_alpha = digamma(a_plus_r);
        if constexpr (!is_constant_all<T_r>::value) {
          partials<0>(ops_partials)[i]
              += digamma(r_plus_n + 1)
                 + (digamma_r_alpha - digamma_n_r_alpha_beta)
                 + (dF[2] + dF[4]) / F - digamma(r_dbl);
        }
        if constexpr (!is_constant_all<T_alpha>::value) {
          partials<1>(ops_partials)[i] += digamma_r_alpha
                                          - digamma_n_r_alpha_beta + dF[4] / F
                                          - digamma(alpha_dbl);
        }
      }

      if constexpr (!is_constant<T_alpha>::value
                    || !is_constant<T_beta>::value) {
        T_partials_return digamma_alpha_beta = digamma(alpha_dbl + beta_dbl);
        if constexpr (!is_constant<T_alpha>::value) {
          partials<1>(ops_partials)[i] += digamma_alpha_beta;
        }
        if constexpr (!is_constant<T_beta>::value) {
          partials<2>(ops_partials)[i]
              += digamma(b_plus_n + 1) - digamma_n_r_alpha_beta
                 + (dF[1] + dF[4]) / F
                 - (digamma(beta_dbl) - digamma_alpha_beta);
        }
      }
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
