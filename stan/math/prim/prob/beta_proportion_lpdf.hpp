#ifndef STAN_MATH_PRIM_PROB_BETA_PROPORTION_LPDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_PROPORTION_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
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
 * @tparam T_y type of scalar outcome
 * @tparam T_loc type of prior location
 * @tparam T_prec type of prior precision
 *
 * @param y (Sequence of) scalar(s) between zero and one
 * @param mu (Sequence of) location parameter(s)
 * @param kappa (Sequence of) precision parameter(s)
 * @return The log of the product of densities.
 */
template <bool propto, typename T_y, typename T_loc, typename T_prec,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_prec>* = nullptr>
return_type_t<T_y, T_loc, T_prec> beta_proportion_lpdf(const T_y& y,
                                                       const T_loc& mu,
                                                       const T_prec& kappa) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_prec>;
  using std::log;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_kappa_ref = ref_type_if_t<!is_constant<T_prec>::value, T_prec>;
  static const char* function = "beta_proportion_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Precision parameter", kappa);
  if (size_zero(y, mu, kappa)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_kappa_ref kappa_ref = kappa;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) kappa_val = to_ref(as_value_column_array_or_scalar(kappa_ref));

  check_positive(function, "Location parameter", mu_val);
  check_less(function, "Location parameter", mu_val, 1.0);
  check_positive_finite(function, "Precision parameter", kappa_val);
  check_bounded(function, "Random variable", value_of(y_val), 0, 1);

  if (!include_summand<propto, T_y, T_loc, T_prec>::value) {
    return 0;
  }

  const auto& log_y
      = to_ref_if<!is_constant_all<T_loc, T_prec>::value>(log(y_val));
  const auto& log1m_y
      = to_ref_if<!is_constant_all<T_loc, T_prec>::value>(log1m(y_val));
  const auto& mukappa = to_ref(mu_val * kappa_val);

  size_t N = max_size(y, mu, kappa);
  T_partials_return logp(0);
  if (include_summand<propto, T_prec>::value) {
    logp += sum(lgamma(kappa_val)) * N / math::size(kappa);
  }
  if (include_summand<propto, T_loc, T_prec>::value) {
    logp -= sum(lgamma(mukappa) + lgamma(kappa_val - mukappa)) * N
            / max_size(mu, kappa_val);
  }
  logp += sum((mukappa - 1) * log_y + (kappa_val - mukappa - 1) * log1m_y);

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, kappa_ref);
  if (!is_constant_all<T_y>::value) {
    edge<0>(ops_partials).partials_
        = (mukappa - 1) / y_val + (kappa_val - mukappa - 1) / (y_val - 1);
  }
  if (!is_constant_all<T_loc, T_prec>::value) {
    auto digamma_mukappa
        = to_ref_if<(!is_constant_all<T_loc>::value
                     && !is_constant_all<T_prec>::value)>(digamma(mukappa));
    auto digamma_kappa_mukappa = to_ref_if<(
        !is_constant_all<T_loc>::value && !is_constant_all<T_prec>::value)>(
        digamma(kappa_val - mukappa));
    if (!is_constant_all<T_loc>::value) {
      edge<1>(ops_partials).partials_
          = kappa_val
            * (digamma_kappa_mukappa - digamma_mukappa + log_y - log1m_y);
    }
    if (!is_constant_all<T_prec>::value) {
      edge<2>(ops_partials).partials_
          = digamma(kappa_val) + mu_val * (log_y - digamma_mukappa)
            + (1 - mu_val) * (log1m_y - digamma_kappa_mukappa);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_prec>
inline return_type_t<T_y, T_loc, T_prec> beta_proportion_lpdf(
    const T_y& y, const T_loc& mu, const T_prec& kappa) {
  return beta_proportion_lpdf<false>(y, mu, kappa);
}

}  // namespace math
}  // namespace stan
#endif
