#ifndef STAN_MATH_PRIM_PROB_BETA_PROPORTION_LCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_PROPORTION_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
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
 * Returns the beta log cumulative distribution function
 * for specified probability, location, and precision parameters:
 * beta_proportion_lcdf(y | mu, kappa) = beta_lcdf(y | mu * kappa, (1 -
 * mu) * kappa).  Any arguments other than scalars must be containers of
 * the same size.  With non-scalar arguments, the return is the sum of
 * the log cdfs with scalars broadcast as necessary.
 *
 * @tparam T_y type of y
 * @tparam T_loc type of location parameter
 * @tparam T_prec type of precision parameter
 * @param y (Sequence of) scalar(s) between zero and one
 * @param mu (Sequence of) location parameter(s)
 * @param kappa (Sequence of) precision parameter(s)
 * @return log probability or sum of log of probabilities
 * @throw std::domain_error if mu is outside of (0, 1)
 * @throw std::domain_error if kappa is nonpositive
 * @throw std::domain_error if y is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_loc, typename T_prec>
return_type_t<T_y, T_loc, T_prec> beta_proportion_lcdf(const T_y& y,
                                                       const T_loc& mu,
                                                       const T_prec& kappa) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_prec>;
  using std::exp;
  using std::log;
  using std::pow;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_kappa_ref = ref_type_t<T_prec>;
  static const char* function = "beta_proportion_lcdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Precision parameter", kappa);
  if (size_zero(y, mu, kappa)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_kappa_ref kappa_ref = kappa;
  check_positive(function, "Location parameter", value_of(mu_ref));
  check_less(function, "Location parameter", value_of(mu_ref), 1.0);
  check_positive_finite(function, "Precision parameter", value_of(kappa_ref));
  check_bounded(function, "Random variable", value_of(y_ref), 0.0, 1.0);

  T_partials_return cdf_log(0.0);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, kappa_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_kappa_ref> kappa_vec(kappa_ref);
  size_t size_kappa = stan::math::size(kappa);
  size_t size_mu_kappa = max_size(mu, kappa);
  size_t N = max_size(y, mu, kappa);

  VectorBuilder<!is_constant_all<T_loc, T_prec>::value, T_partials_return,
                T_loc, T_prec>
      digamma_mukappa(size_mu_kappa);
  VectorBuilder<!is_constant_all<T_loc, T_prec>::value, T_partials_return,
                T_loc, T_prec>
      digamma_kappa_mukappa(size_mu_kappa);
  VectorBuilder<!is_constant_all<T_loc, T_prec>::value, T_partials_return,
                T_prec>
      digamma_kappa(size_kappa);

  if (!is_constant_all<T_loc, T_prec>::value) {
    for (size_t i = 0; i < size_mu_kappa; i++) {
      const T_partials_return kappa_dbl = kappa_vec.val(i);
      const T_partials_return mukappa_dbl = mu_vec.val(i) * kappa_dbl;
      digamma_mukappa[i] = digamma(mukappa_dbl);
      digamma_kappa_mukappa[i] = digamma(kappa_dbl - mukappa_dbl);
    }
    for (size_t i = 0; i < size_kappa; i++) {
      digamma_kappa[i] = digamma(kappa_vec.val(i));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return mu_dbl = mu_vec.val(n);
    const T_partials_return kappa_dbl = kappa_vec.val(n);
    const T_partials_return mukappa_dbl = mu_dbl * kappa_dbl;
    const T_partials_return kappa_mukappa_dbl = kappa_dbl - mukappa_dbl;
    const T_partials_return betafunc_dbl = beta(mukappa_dbl, kappa_mukappa_dbl);
    const T_partials_return Pn
        = inc_beta(mukappa_dbl, kappa_mukappa_dbl, y_dbl);

    cdf_log += log(Pn);

    const T_partials_return inv_Pn
        = is_constant_all<T_y, T_loc, T_prec>::value ? 0 : inv(Pn);

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n] += pow(1 - y_dbl, kappa_mukappa_dbl - 1)
                                      * pow(y_dbl, mukappa_dbl - 1) * inv_Pn
                                      / betafunc_dbl;
    }

    T_partials_return g1 = 0;
    T_partials_return g2 = 0;

    if (!is_constant_all<T_loc, T_prec>::value) {
      grad_reg_inc_beta(g1, g2, mukappa_dbl, kappa_mukappa_dbl, y_dbl,
                        digamma_mukappa[n], digamma_kappa_mukappa[n],
                        digamma_kappa[n], betafunc_dbl);
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials)[n] += kappa_dbl * (g1 - g2) * inv_Pn;
    }
    if (!is_constant_all<T_prec>::value) {
      partials<2>(ops_partials)[n]
          += (g1 * mu_dbl + g2 * (1 - mu_dbl)) * inv_Pn;
    }
  }

  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
