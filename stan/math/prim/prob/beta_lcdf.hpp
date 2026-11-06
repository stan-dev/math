#ifndef STAN_MATH_PRIM_PROB_BETA_LCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

#include <cmath>

namespace stan {
namespace math {

/**
 * Beta log CDF.
 *
 * @tparam T_y type of y
 * @tparam T_scale_succ type of success parameter
 * @tparam T_scale_fail type of failure parameter
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline return_type_t<T_y, T_scale_succ, T_scale_fail> beta_lcdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta_param) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale_succ>;
  using T_beta_ref = ref_type_t<T_scale_fail>;
  static constexpr const char* function = "beta_lcdf";

  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta_param);
  if (size_zero(y, alpha, beta_param)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta_param;
  check_positive_finite(function, "First shape parameter", alpha_ref);
  check_positive_finite(function, "Second shape parameter", beta_ref);
  check_bounded(function, "Random variable", value_of(y_ref), 0, 1);

  T_partials_return cdf_log(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);

  const size_t N = max_size(y, alpha, beta_param);

  constexpr bool calc_alpha = is_autodiff_v<T_scale_succ>;
  constexpr bool calc_beta = is_autodiff_v<T_scale_fail>;
  for (size_t n = 0; n < N; ++n) {
    const auto y_dbl = y_vec.val(n);
    const auto alpha_dbl = alpha_vec.val(n);
    const auto beta_dbl = beta_vec.val(n);
    const auto betafunc_dbl = beta(alpha_dbl, beta_dbl);
    const auto Pn = inc_beta(alpha_dbl, beta_dbl, y_dbl);

    const auto inv_Pn
        = is_any_autodiff_v<T_y, T_scale_succ, T_scale_fail> ? inv(Pn) : 0;

    cdf_log += log(Pn);

    if constexpr (is_any_autodiff_v<T_y>) {
      partials<0>(ops_partials)[n] += pow(1 - y_dbl, beta_dbl - 1)
                                      * pow(y_dbl, alpha_dbl - 1) * inv_Pn
                                      / betafunc_dbl;
    }

    if constexpr (calc_alpha && calc_beta) {
      const auto [g1, g2] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
          alpha_dbl, beta_dbl, y_dbl, digamma(alpha_dbl), digamma(beta_dbl),
          digamma(alpha_dbl + beta_dbl), 0.0);
      partials<1>(ops_partials)[n] += g1 * inv_Pn;
      partials<2>(ops_partials)[n] += g2 * inv_Pn;
    } else if constexpr (calc_alpha) {
      partials<1>(ops_partials)[n]
          += grad_reg_inc_beta<GRAD_RIB::G1>(alpha_dbl, beta_dbl, y_dbl,
                                             digamma(alpha_dbl), 0.0,
                                             digamma(alpha_dbl + beta_dbl), 0.0)
             * inv_Pn;
    } else if constexpr (calc_beta) {
      partials<2>(ops_partials)[n]
          += grad_reg_inc_beta<GRAD_RIB::G2>(alpha_dbl, beta_dbl, y_dbl, 0.0,
                                             digamma(beta_dbl),
                                             digamma(alpha_dbl + beta_dbl), 0.0)
             * inv_Pn;
    }
  }

  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif  // STAN_MATH_PRIM_PROB_BETA_LCDF_HPP
