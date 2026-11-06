#ifndef STAN_MATH_PRIM_PROB_BETA_CDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the beta cumulative distribution function for the given
 * variate and scale variables.
 *
 * @param y A scalar variate.
 * @param alpha Prior sample size.
 * @param beta Prior sample size.
 * @return The beta cdf evaluated at the specified arguments.
 * @tparam T_y Type of y.
 * @tparam T_scale_succ Type of alpha.
 * @tparam T_scale_fail Type of beta.
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline return_type_t<T_y, T_scale_succ, T_scale_fail> beta_cdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale_succ>;
  using T_beta_ref = ref_type_t<T_scale_fail>;
  static constexpr const char* function = "beta_cdf";
  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta);
  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "First shape parameter", alpha_ref);
  check_positive_finite(function, "Second shape parameter", beta_ref);
  check_bounded(function, "Random variable", value_of(y_ref), 0, 1);

  T_partials_return P(1.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) <= 0) {
      return ops_partials.build(0.0);
    }
  }

  constexpr bool calc_alpha = is_autodiff_v<T_scale_succ>;
  constexpr bool calc_beta = is_autodiff_v<T_scale_fail>;
  for (size_t n = 0; n < N; n++) {
    const auto y_dbl = y_vec.val(n);

    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_dbl >= 1.0) {
      continue;
    }

    const auto alpha_dbl = alpha_vec.val(n);
    const auto beta_dbl = beta_vec.val(n);
    const auto Pn = inc_beta(alpha_dbl, beta_dbl, y_dbl);

    P *= Pn;
    if constexpr (is_any_autodiff_v<T_y, T_scale_succ, T_scale_fail>) {
      const auto inv_Pn = inv(Pn);
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n]
            += inc_beta_ddz(alpha_dbl, beta_dbl, y_dbl) * inv_Pn;
      }

      if constexpr (calc_alpha && calc_beta) {
        const auto [g1, g2] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
            alpha_dbl, beta_dbl, y_dbl, digamma(alpha_dbl), digamma(beta_dbl),
            digamma(alpha_dbl + beta_dbl), 0.0);
        partials<1>(ops_partials)[n] += g1 * inv_Pn;
        partials<2>(ops_partials)[n] += g2 * inv_Pn;
      } else if constexpr (calc_alpha) {
        partials<1>(ops_partials)[n]
            += grad_reg_inc_beta<GRAD_RIB::G1>(
                   alpha_dbl, beta_dbl, y_dbl, digamma(alpha_dbl), 0.0,
                   digamma(alpha_dbl + beta_dbl), 0.0)
               * inv_Pn;
      } else if constexpr (calc_beta) {
        partials<2>(ops_partials)[n]
            += grad_reg_inc_beta<GRAD_RIB::G2>(
                   alpha_dbl, beta_dbl, y_dbl, 0.0, digamma(beta_dbl),
                   digamma(alpha_dbl + beta_dbl), 0.0)
               * inv_Pn;
      }
    }
  }

  if constexpr (is_autodiff_v<T_y>) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_scale_succ>) {
    for (size_t n = 0; n < stan::math::size(alpha); ++n) {
      partials<1>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_scale_fail>) {
    for (size_t n = 0; n < stan::math::size(beta); ++n) {
      partials<2>(ops_partials)[n] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
