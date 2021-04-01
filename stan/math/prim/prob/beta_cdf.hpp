#ifndef STAN_MATH_PRIM_PROB_BETA_CDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddb.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
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
return_type_t<T_y, T_scale_succ, T_scale_fail> beta_cdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale_succ>;
  using T_beta_ref = ref_type_t<T_scale_fail>;
  static const char* function = "beta_cdf";
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
  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(
      y_ref, alpha_ref, beta_ref);
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_alpha = stan::math::size(alpha);
  size_t size_beta = stan::math::size(beta);
  size_t size_alpha_beta = max_size(alpha, beta);
  size_t N = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) <= 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_scale_succ>::value, T_partials_return,
                T_scale_succ>
      digamma_alpha(size_alpha);
  if (!is_constant_all<T_scale_succ>::value) {
    for (size_t n = 0; n < size_alpha; n++) {
      digamma_alpha[n] = digamma(alpha_vec.val(n));
    }
  }

  VectorBuilder<!is_constant_all<T_scale_fail>::value, T_partials_return,
                T_scale_fail>
      digamma_beta(size_beta);
  if (!is_constant_all<T_scale_fail>::value) {
    for (size_t n = 0; n < size_beta; n++) {
      digamma_beta[n] = digamma(beta_vec.val(n));
    }
  }

  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_succ, T_scale_fail>
      digamma_sum(size_alpha_beta);
  if (!is_constant_all<T_scale_succ, T_scale_fail>::value) {
    for (size_t n = 0; n < size_alpha_beta; n++) {
      digamma_sum[n] = digamma(alpha_vec.val(n) + beta_vec.val(n));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);

    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_dbl >= 1.0) {
      continue;
    }

    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return beta_dbl = beta_vec.val(n);
    const T_partials_return Pn = inc_beta(alpha_dbl, beta_dbl, y_dbl);
    const T_partials_return inv_Pn
        = is_constant_all<T_y, T_scale_succ, T_scale_fail>::value ? 0 : inv(Pn);

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += inc_beta_ddz(alpha_dbl, beta_dbl, y_dbl) * inv_Pn;
    }

    if (!is_constant_all<T_scale_succ>::value) {
      ops_partials.edge2_.partials_[n]
          += inc_beta_dda(alpha_dbl, beta_dbl, y_dbl, digamma_alpha[n],
                          digamma_sum[n])
             * inv_Pn;
    }
    if (!is_constant_all<T_scale_fail>::value) {
      ops_partials.edge3_.partials_[n]
          += inc_beta_ddb(alpha_dbl, beta_dbl, y_dbl, digamma_beta[n],
                          digamma_sum[n])
             * inv_Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale_succ>::value) {
    for (size_t n = 0; n < stan::math::size(alpha); ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale_fail>::value) {
    for (size_t n = 0; n < stan::math::size(beta); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
