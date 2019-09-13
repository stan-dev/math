#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The log of the beta density for specified y, location, and
 * precision: beta_proportion_lpdf(y | mu, kappa) = beta_lpdf(y | mu *
 * kappa, (1 - mu) * kappa).  Any arguments other than scalars must be
 * containers of the same size.  With non-scalar arguments, the return
 * is the sum of the log pdfs with scalars broadcast as necessary.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/mu/kappa triple.
 *
 * Prior location, mu, must be contained in (0, 1).  Prior precision
 * must be positive.
 *
 * @param y (Sequence of) scalar(s) between zero and one
 * @param mu (Sequence of) location parameter(s)
 * @param kappa (Sequence of) precision parameter(s)
 * @return The log of the product of densities.
 * @tparam T_y Type of scalar outcome.
 * @tparam T_loc Type of prior location.
 * @tparam T_prec Type of prior precision.
 */
template <bool propto, typename T_y, typename T_loc, typename T_prec>
inline auto beta_proportion_lpdf(T_y&& y, T_loc&& mu,
                                 T_prec&& kappa) {
  using T_partials = partials_return_t<T_y, T_loc, T_prec>;
  T_partials logp(0);
  using std::log;

  static const char* function = "beta_proportion_lpdf";
  check_positive(function, "Location parameter", mu);
  check_less_or_equal(function, "Location parameter", mu, 1.0);
  check_positive_finite(function, "Precision parameter", kappa);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_less_or_equal(function, "Random variable", y, 1.0);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Precision parameter", kappa);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_loc> mu_vec(mu);
  const scalar_seq_view<T_prec> kappa_vec(kappa);
  const size_t N = max_size(y, mu, kappa);
  const size_t N_mukappa = max_size(mu, kappa);
  operands_and_partials<T_y, T_loc, T_prec> ops_partials(y, mu, kappa);
  if (!include_summand<propto, T_y, T_loc, T_prec>::value) {
    return ops_partials.build(logp);
  } else if (size_zero(y, mu, kappa)) {
    return ops_partials.build(logp);
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    if (y_dbl < 0 || y_dbl > 1) {
      return ops_partials.build(T_partials(LOG_ZERO));
    }
  }
  VectorBuilder<include_summand<propto, T_y, T_loc, T_prec>::value, T_partials,
                T_y>
      log_y(length(y));
  VectorBuilder<include_summand<propto, T_y, T_loc, T_prec>::value, T_partials,
                T_y>
      log1m_y(length(y));

  for (size_t n = 0; n < length(y); n++) {
    if (include_summand<propto, T_y, T_loc, T_prec>::value) {
      log_y[n] = log(value_of(y_vec[n]));
      log1m_y[n] = log1m(value_of(y_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_loc, T_prec>::value, T_partials,
                T_loc, T_prec>
      lgamma_mukappa(N_mukappa);
  VectorBuilder<include_summand<propto, T_loc, T_prec>::value, T_partials,
                T_loc, T_prec>
      lgamma_kappa_mukappa(N_mukappa);
  VectorBuilder<!is_constant_all<T_loc, T_prec>::value, T_partials, T_loc,
                T_prec>
      digamma_mukappa(N_mukappa);
  VectorBuilder<!is_constant_all<T_loc, T_prec>::value, T_partials, T_loc,
                T_prec>
      digamma_kappa_mukappa(N_mukappa);

  for (size_t n = 0; n < N_mukappa; n++) {
    const T_partials mukappa_dbl = value_of(mu_vec[n]) * value_of(kappa_vec[n]);
    const T_partials kappa_mukappa_dbl = value_of(kappa_vec[n]) - mukappa_dbl;

    if (include_summand<propto, T_loc, T_prec>::value) {
      lgamma_mukappa[n] = lgamma(mukappa_dbl);
      lgamma_kappa_mukappa[n] = lgamma(kappa_mukappa_dbl);
    }

    if (!is_constant_all<T_loc, T_prec>::value) {
      digamma_mukappa[n] = digamma(mukappa_dbl);
      digamma_kappa_mukappa[n] = digamma(kappa_mukappa_dbl);
    }
  }

  VectorBuilder<include_summand<propto, T_prec>::value, T_partials, T_prec>
      lgamma_kappa(length(kappa));
  VectorBuilder<!is_constant_all<T_prec>::value, T_partials, T_prec>
      digamma_kappa(length(kappa));

  for (size_t n = 0; n < length(kappa); n++) {
    if (include_summand<propto, T_prec>::value) {
      lgamma_kappa[n] = lgamma(value_of(kappa_vec[n]));
    }

    if (!is_constant_all<T_prec>::value) {
      digamma_kappa[n] = digamma(value_of(kappa_vec[n]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    const T_partials mu_dbl = value_of(mu_vec[n]);
    const T_partials kappa_dbl = value_of(kappa_vec[n]);

    if (include_summand<propto, T_prec>::value) {
      logp += lgamma_kappa[n];
    }
    if (include_summand<propto, T_loc, T_prec>::value) {
      logp -= lgamma_mukappa[n] + lgamma_kappa_mukappa[n];
    }
    if (include_summand<propto, T_y, T_loc, T_prec>::value) {
      const T_partials mukappa_dbl = mu_dbl * kappa_dbl;
      logp += (mukappa_dbl - 1) * log_y[n]
              + (kappa_dbl - mukappa_dbl - 1) * log1m_y[n];
    }

    if (!is_constant_all<T_y>::value) {
      const T_partials mukappa_dbl = mu_dbl * kappa_dbl;
      ops_partials.edge1_.partials_[n]
          += (mukappa_dbl - 1) / y_dbl
             + (kappa_dbl - mukappa_dbl - 1) / (y_dbl - 1);
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += kappa_dbl
             * (digamma_kappa_mukappa[n] - digamma_mukappa[n] + log_y[n]
                - log1m_y[n]);
    }
    if (!is_constant_all<T_prec>::value) {
      ops_partials.edge3_.partials_[n]
          += digamma_kappa[n] + mu_dbl * (log_y[n] - digamma_mukappa[n])
             + (1 - mu_dbl) * (log1m_y[n] - digamma_kappa_mukappa[n]);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_prec>
inline auto beta_proportion_lpdf(T_y&& y, T_loc&& mu,
                                 T_prec&& kappa) {
  return beta_proportion_lpdf<false>(y, mu, kappa);
}

}  // namespace math
}  // namespace stan
#endif
