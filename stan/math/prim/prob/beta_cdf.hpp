#ifndef STAN_MATH_PRIM_PROB_BETA_CDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddb.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
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

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  static const char* function = "beta_cdf";

  T_partials_return P(1.0);

  check_positive_finite(function, "First shape parameter", alpha);
  check_positive_finite(function, "Second shape parameter", beta);
  check_not_nan(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta);
  check_nonnegative(function, "Random variable", y);
  check_less_or_equal(function, "Random variable", y, 1);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_scale_succ> alpha_vec(alpha);
  scalar_seq_view<T_scale_fail> beta_vec(beta);
  size_t N = max_size(y, alpha, beta);

  operands_and_partials<T_y, T_scale_succ, T_scale_fail> ops_partials(y, alpha,
                                                                      beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(y); i++) {
    if (value_of(y_vec[i]) <= 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_succ, T_scale_fail>
      digamma_alpha_vec(max_size(alpha, beta));

  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_succ, T_scale_fail>
      digamma_beta_vec(max_size(alpha, beta));

  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_succ, T_scale_fail>
      digamma_sum_vec(max_size(alpha, beta));

  if (!is_constant_all<T_scale_succ, T_scale_fail>::value) {
    for (size_t n = 0; n < N; n++) {
      const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
      const T_partials_return beta_dbl = value_of(beta_vec[n]);

      digamma_alpha_vec[n] = digamma(alpha_dbl);
      digamma_beta_vec[n] = digamma(beta_dbl);
      digamma_sum_vec[n] = digamma(alpha_dbl + beta_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) >= 1.0) {
      continue;
    }

    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return beta_dbl = value_of(beta_vec[n]);

    const T_partials_return Pn = inc_beta(alpha_dbl, beta_dbl, y_dbl);

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += inc_beta_ddz(alpha_dbl, beta_dbl, y_dbl) / Pn;
    }

    if (!is_constant_all<T_scale_succ>::value) {
      ops_partials.edge2_.partials_[n]
          += inc_beta_dda(alpha_dbl, beta_dbl, y_dbl, digamma_alpha_vec[n],
                          digamma_sum_vec[n])
             / Pn;
    }
    if (!is_constant_all<T_scale_fail>::value) {
      ops_partials.edge3_.partials_[n]
          += inc_beta_ddb(alpha_dbl, beta_dbl, y_dbl, digamma_beta_vec[n],
                          digamma_sum_vec[n])
             / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale_succ>::value) {
    for (size_t n = 0; n < size(alpha); ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale_fail>::value) {
    for (size_t n = 0; n < size(beta); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
