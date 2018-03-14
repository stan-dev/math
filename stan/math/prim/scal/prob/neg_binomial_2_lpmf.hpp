#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LPMF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LPMF_HPP

#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/prob/poisson_lpmf.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/random/negative_binomial_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace stan {
namespace math {

// NegBinomial(n|mu, phi)  [mu >= 0; phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_location, typename T_precision>
typename return_type<T_location, T_precision>::type neg_binomial_2_lpmf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  typedef
      typename stan::partials_return_type<T_n, T_location, T_precision>::type
          T_partials_return;

  static const char* function = "neg_binomial_2_lpmf";

  if (size_zero(n, mu, phi))
    return 0.0;

  T_partials_return logp(0.0);
  check_nonnegative(function, "Failures variable", n);
  check_positive_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Precision parameter", phi);
  check_consistent_sizes(function, "Failures variable", n, "Location parameter",
                         mu, "Precision parameter", phi);

  if (!include_summand<propto, T_location, T_precision>::value)
    return 0.0;

  using std::log;
  using std::log;

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_location> mu_vec(mu);
  scalar_seq_view<T_precision> phi_vec(phi);
  size_t size = max_size(n, mu, phi);

  operands_and_partials<T_location, T_precision> ops_partials(mu, phi);

  size_t len_ep = max_size(mu, phi);
  size_t len_np = max_size(n, phi);

  VectorBuilder<true, T_partials_return, T_location> mu__(length(mu));
  #pragma omp parallel for if (length(mu) > 3 * omp_get_max_threads()) \
    default(none) shared(mu__, mu_vec, mu)
  for (size_t i = 0; i < length(mu); ++i)
    mu__[i] = value_of(mu_vec[i]);

  VectorBuilder<true, T_partials_return, T_precision> phi__(length(phi));
  VectorBuilder<true, T_partials_return, T_precision> log_phi(length(phi));
  #pragma omp parallel for if (length(phi) > 3 * omp_get_max_threads()) \
    default(none) shared(phi__, phi_vec, log_phi, phi)
  for (size_t i = 0; i < length(phi); ++i) {
    phi__[i] = value_of(phi_vec[i]);
    log_phi[i] = log(phi__[i]);
  }

  VectorBuilder<true, T_partials_return, T_location, T_precision>
      log_mu_plus_phi(len_ep);
  #pragma omp parallel for if (len_ep > 3 * omp_get_max_threads()) \
    default(none) shared(log_mu_plus_phi, mu__, phi__, len_ep)
  for (size_t i = 0; i < len_ep; ++i)
    log_mu_plus_phi[i] = log(mu__[i] + phi__[i]);

  VectorBuilder<true, T_partials_return, T_n, T_precision> n_plus_phi(len_np);
  #pragma omp parallel for if (len_np > 3 * omp_get_max_threads()) \
    default(none) shared(n_plus_phi, n_vec, phi__, len_np)
  for (size_t i = 0; i < len_np; ++i)
    n_plus_phi[i] = n_vec[i] + phi__[i];

#ifndef STAN_MATH_FWD_CORE_HPP
  #pragma omp parallel for if (size > 3 * omp_get_max_threads()) \
    reduction(+ : logp) default(none) \
    shared(n_vec, phi__, n_plus_phi, log_mu_plus_phi, mu__, ops_partials, \
           log_phi, size)
#endif
  for (size_t i = 0; i < size; i++) {
    if (include_summand<propto>::value)
      logp -= lgamma(n_vec[i] + 1.0);
    if (include_summand<propto, T_precision>::value)
      logp += multiply_log(phi__[i], phi__[i]) - lgamma(phi__[i]);
    if (include_summand<propto, T_location, T_precision>::value)
      logp -= (n_plus_phi[i]) * log_mu_plus_phi[i];
    if (include_summand<propto, T_location>::value)
      logp += multiply_log(n_vec[i], mu__[i]);
    if (include_summand<propto, T_precision>::value)
      logp += lgamma(n_plus_phi[i]);

    // if phi is large we probably overflow, defer to Poisson:
    if (phi__[i] > 1e5) {
      logp = poisson_lpmf(n_vec[i], mu__[i]); // should this be += ?
    }

    if (!is_constant_struct<T_location>::value)
      ops_partials.edge1_.partials_[i]
          += n_vec[i] / mu__[i] - (n_vec[i] + phi__[i]) / (mu__[i] + phi__[i]);
    if (!is_constant_struct<T_precision>::value)
      ops_partials.edge2_.partials_[i]
          += 1.0 - n_plus_phi[i] / (mu__[i] + phi__[i]) + log_phi[i]
             - log_mu_plus_phi[i] - digamma(phi__[i]) + digamma(n_plus_phi[i]);
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_location, typename T_precision>
inline typename return_type<T_location, T_precision>::type neg_binomial_2_lpmf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  return neg_binomial_2_lpmf<false>(n, mu, phi);
}

}  // namespace math
}  // namespace stan
#endif
