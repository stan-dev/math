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
 * @param precision precision for `grad_F32`, default \f$10^{-8}\f$
 * @param max_steps max iteration allowed for `grad_F32`, default \f$10^{-8}\f$
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if r, alpha, or beta fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_r, typename T_alpha, typename T_beta>
inline return_type_t<T_r, T_alpha, T_beta> beta_neg_binomial_lccdf(
    const T_n& n, const T_r& r, const T_alpha& alpha, const T_beta& beta,
    const double precision = 1e-8, const int max_steps = 1e6) {
  static constexpr const char* function = "beta_neg_binomial_lccdf";
  check_consistent_sizes(
      function, "Failures variable", n, "Number of successes parameter", r,
      "Prior success parameter", alpha, "Prior failure parameter", beta);
  if (size_zero(n, r, alpha, beta)) {
    return 0;
  }

  using T_r_ref = ref_type_t<T_r>;
  T_r_ref r_ref = r;
  using T_alpha_ref = ref_type_t<T_alpha>;
  T_alpha_ref alpha_ref = alpha;
  using T_beta_ref = ref_type_t<T_beta>;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Number of successes parameter", r_ref);
  check_positive_finite(function, "Prior success parameter", alpha_ref);
  check_positive_finite(function, "Prior failure parameter", beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_r_ref> r_vec(r_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  int size_n = stan::math::size(n);
  size_t max_size_seq_view = max_size(n, r, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (int i = 0; i < size_n; i++) {
    if (n_vec.val(i) < 0) {
      return 0.0;
    }
  }

  using T_partials_return = partials_return_t<T_n, T_r, T_alpha, T_beta>;
  T_partials_return log_ccdf(0.0);
  auto ops_partials = make_partials_propagator(r_ref, alpha_ref, beta_ref);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit return for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(negative_infinity());
    }
    auto n_dbl = n_vec.val(i);
    auto r_dbl = r_vec.val(i);
    auto alpha_dbl = alpha_vec.val(i);
    auto beta_dbl = beta_vec.val(i);
    auto b_plus_n = beta_dbl + n_dbl;
    auto r_plus_n = r_dbl + n_dbl;
    auto a_plus_r = alpha_dbl + r_dbl;
    using a_t = return_type_t<decltype(b_plus_n), decltype(r_plus_n)>;
    using b_t = return_type_t<decltype(n_dbl), decltype(a_plus_r),
                              decltype(b_plus_n)>;
    auto F = hypergeometric_3F2(
        std::initializer_list<a_t>{1.0, b_plus_n + 1.0, r_plus_n + 1.0},
        std::initializer_list<b_t>{n_dbl + 2.0, a_plus_r + b_plus_n + 1.0},
        1.0);
    auto C = lgamma(r_plus_n + 1.0) + lbeta(a_plus_r, b_plus_n + 1.0)
             - lgamma(r_dbl) - lbeta(alpha_dbl, beta_dbl) - lgamma(n_dbl + 2);
    log_ccdf += C + stan::math::log(F);

    if constexpr (!is_constant_all<T_r, T_alpha, T_beta>::value) {
      auto digamma_n_r_alpha_beta = digamma(a_plus_r + b_plus_n + 1.0);
      T_partials_return dF[6];
      grad_F32<false, !is_constant<T_beta>::value, !is_constant_all<T_r>::value,
               false, true, false>(dF, 1.0, b_plus_n + 1.0, r_plus_n + 1.0,
                                   n_dbl + 2.0, a_plus_r + b_plus_n + 1.0, 1.0,
                                   precision, max_steps);

      if constexpr (!is_constant<T_r>::value || !is_constant<T_alpha>::value) {
        auto digamma_r_alpha = digamma(a_plus_r);
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
        auto digamma_alpha_beta = digamma(alpha_dbl + beta_dbl);
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

  return ops_partials.build(log_ccdf);
}

}  // namespace math
}  // namespace stan
#endif
