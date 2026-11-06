#ifndef STAN_MATH_PRIM_PROB_BETA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the beta log complementary cumulative distribution function
 * for the given probability, success, and failure parameters.  Any
 * arguments other than scalars must be containers of the same size.
 * With non-scalar arguments, the return is the sum of the log ccdfs
 * with scalars broadcast as necessary.
 *
 * @tparam T_y type of y
 * @tparam T_scale_succ type of success parameter
 * @tparam T_scale_fail type of failure parameter
 * @param y (Sequence of) scalar(s) between zero and one
 * @param alpha (Sequence of) success parameter(s)
 * @param beta_param (Sequence of) failure parameter(s)
 * @return log probability or sum of log of probabilities
 * @throw std::domain_error if alpha or beta is negative
 * @throw std::domain_error if y is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline return_type_t<T_y, T_scale_succ, T_scale_fail> beta_lccdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta_param) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using std::exp;
  using std::log;
  using std::pow;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale_succ>;
  using T_beta_ref = ref_type_t<T_scale_fail>;
  static constexpr const char* function = "beta_lccdf";
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

  T_partials_return ccdf_log(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(y, alpha, beta_param);

  constexpr bool calc_alpha = is_autodiff_v<T_scale_succ>;
  constexpr bool calc_beta = is_autodiff_v<T_scale_fail>;
  for (size_t n = 0; n < N; n++) {
    const auto y_dbl = y_vec.val(n);

    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_dbl >= 1.0) {
      return ops_partials.build(negative_infinity());
    }

    const auto alpha_dbl = alpha_vec.val(n);
    const auto beta_dbl = beta_vec.val(n);
    // The complement 1 - I_y(a, b), evaluated from the smaller of y and
    // 1 - y. Forming it as a difference loses every digit once I_y(a, b)
    // rounds to 1, which happens for every complement below eps.
    const auto Pn = internal::inc_beta_complement(alpha_dbl, beta_dbl, y_dbl,
                                                  1.0 - y_dbl);
    const auto inv_Pn
        = is_constant_all<T_y, T_scale_succ, T_scale_fail>::value ? 0 : inv(Pn);

    ccdf_log += log(Pn);

    if constexpr (is_autodiff_v<T_y>) {
      // inc_beta_ddz is the same density, evaluated without the product of
      // two powers over a beta function, each of which can underflow on
      // its own.
      partials<0>(ops_partials)[n]
          -= inc_beta_ddz(alpha_dbl, beta_dbl, y_dbl) * inv_Pn;
    }

    // The shape derivatives from the same argument as the value. For
    // y <= 1/2 they are those of I_y(a, b), which are those of -Pn; on the
    // reflected arguments they are those of Pn itself, with beta the first
    // shape and alpha the second.
    [[maybe_unused]] const bool reflect = y_dbl > 0.5;
    if constexpr (calc_alpha && calc_beta) {
      if (reflect) {
        const auto [g1, g2] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
            beta_dbl, alpha_dbl, 1.0 - y_dbl, digamma(beta_dbl),
            digamma(alpha_dbl), digamma(alpha_dbl + beta_dbl), 0.0);
        partials<1>(ops_partials)[n] += g2 * inv_Pn;
        partials<2>(ops_partials)[n] += g1 * inv_Pn;
      } else {
        const auto [g1, g2] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
            alpha_dbl, beta_dbl, y_dbl, digamma(alpha_dbl), digamma(beta_dbl),
            digamma(alpha_dbl + beta_dbl), 0.0);
        partials<1>(ops_partials)[n] -= g1 * inv_Pn;
        partials<2>(ops_partials)[n] -= g2 * inv_Pn;
      }
    } else if constexpr (calc_alpha) {
      const auto digamma_sum = digamma(alpha_dbl + beta_dbl);
      partials<1>(ops_partials)[n]
          += (reflect ? grad_reg_inc_beta<GRAD_RIB::G2>(
                            beta_dbl, alpha_dbl, 1.0 - y_dbl, 0.0,
                            digamma(alpha_dbl), digamma_sum, 0.0)
                      : -grad_reg_inc_beta<GRAD_RIB::G1>(
                            alpha_dbl, beta_dbl, y_dbl, digamma(alpha_dbl),
                            0.0, digamma_sum, 0.0))
             * inv_Pn;
    } else if constexpr (calc_beta) {
      const auto digamma_sum = digamma(alpha_dbl + beta_dbl);
      partials<2>(ops_partials)[n]
          += (reflect ? grad_reg_inc_beta<GRAD_RIB::G1>(
                            beta_dbl, alpha_dbl, 1.0 - y_dbl,
                            digamma(beta_dbl), 0.0, digamma_sum, 0.0)
                      : -grad_reg_inc_beta<GRAD_RIB::G2>(
                            alpha_dbl, beta_dbl, y_dbl, 0.0, digamma(beta_dbl),
                            digamma_sum, 0.0))
             * inv_Pn;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
