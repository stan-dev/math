#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
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
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/random/negative_binomial_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#ifndef OMP_TRIGGER
#define OMP_TRIGGER 3
#endif
#endif

namespace stan {
namespace math {

// NegBinomial(n|eta, phi)  [phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_log_location,
          typename T_precision>
typename return_type<T_log_location, T_precision>::type neg_binomial_2_log_lpmf(
    const T_n& n, const T_log_location& eta, const T_precision& phi) {
  typedef
      typename stan::partials_return_type<T_n, T_log_location,
                                          T_precision>::type T_partials_return;

  static const char* function = "neg_binomial_2_log_lpmf";

  if (size_zero(n, eta, phi))
    return 0.0;

  T_partials_return logp(0.0);
  check_nonnegative(function, "Failures variable", n);
  check_finite(function, "Log location parameter", eta);
  check_positive_finite(function, "Precision parameter", phi);
  check_consistent_sizes(function, "Failures variable", n,
                         "Log location parameter", eta, "Precision parameter",
                         phi);

  if (!include_summand<propto, T_log_location, T_precision>::value)
    return 0.0;

  using std::exp;
  using std::log;

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_log_location> eta_vec(eta);
  scalar_seq_view<T_precision> phi_vec(phi);
  size_t size = max_size(n, eta, phi);

  operands_and_partials<T_log_location, T_precision> ops_partials(eta, phi);

  size_t len_ep = max_size(eta, phi);
  size_t len_np = max_size(n, phi);

  VectorBuilder<true, T_partials_return, T_log_location> eta__(length(eta));
#pragma omp parallel for if (length(eta) > OMP_TRIGGER * omp_get_max_threads()) default( \
    none) shared(eta__, eta_vec, eta)
  for (size_t i = 0; i < length(eta); ++i)
    eta__[i] = value_of(eta_vec[i]);

  VectorBuilder<true, T_partials_return, T_precision> phi__(length(phi));
#pragma omp parallel for if (length(phi) > OMP_TRIGGER * omp_get_max_threads()) default( \
    none) shared(phi__, phi_vec, phi)
  for (size_t i = 0; i < length(phi); ++i)
    phi__[i] = value_of(phi_vec[i]);

  VectorBuilder<true, T_partials_return, T_precision> log_phi(length(phi));
#pragma omp parallel for if (length(phi) > OMP_TRIGGER * omp_get_max_threads()) default( \
    none) shared(log_phi, phi__, phi)
  for (size_t i = 0; i < length(phi); ++i)
    log_phi[i] = log(phi__[i]);

  VectorBuilder<true, T_partials_return, T_log_location, T_precision>
      logsumexp_eta_logphi(len_ep);
#pragma omp parallel for if (len_ep > OMP_TRIGGER * omp_get_max_threads()) default(none) \
    shared(logsumexp_eta_logphi, eta__, log_phi, len_ep)
  for (size_t i = 0; i < len_ep; ++i)
    logsumexp_eta_logphi[i] = log_sum_exp(eta__[i], log_phi[i]);

  VectorBuilder<true, T_partials_return, T_n, T_precision> n_plus_phi(len_np);
#pragma omp parallel for if (len_np > OMP_TRIGGER * omp_get_max_threads()) default(none) \
    shared(n_plus_phi, n_vec, phi__, len_np)
  for (size_t i = 0; i < len_np; ++i)
    n_plus_phi[i] = n_vec[i] + phi__[i];

#ifndef STAN_MATH_FWD_CORE_HPP
#pragma omp parallel for if (size > OMP_TRIGGER * omp_get_max_threads()) \
    reduction(+ : logp) default(none) \
    shared(n_vec, phi__, n_plus_phi, ops_partials, logsumexp_eta_logphi, \
           eta__, log_phi, size)
#endif
  for (size_t i = 0; i < size; i++) {
    if (include_summand<propto>::value)
      logp -= lgamma(n_vec[i] + 1.0);
    if (include_summand<propto, T_precision>::value)
      logp += multiply_log(phi__[i], phi__[i]) - lgamma(phi__[i]);
    if (include_summand<propto, T_log_location, T_precision>::value)
      logp -= (n_plus_phi[i]) * logsumexp_eta_logphi[i];
    if (include_summand<propto, T_log_location>::value)
      logp += n_vec[i] * eta__[i];
    if (include_summand<propto, T_precision>::value)
      logp += lgamma(n_plus_phi[i]);

    if (!is_constant_struct<T_log_location>::value)
      ops_partials.edge1_.partials_[i]
          += n_vec[i] - n_plus_phi[i] / (phi__[i] / exp(eta__[i]) + 1.0);
    if (!is_constant_struct<T_precision>::value)
      ops_partials.edge2_.partials_[i]
          += 1.0 - n_plus_phi[i] / (exp(eta__[i]) + phi__[i]) + log_phi[i]
             - logsumexp_eta_logphi[i] - digamma(phi__[i])
             + digamma(n_plus_phi[i]);
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_log_location, typename T_precision>
inline typename return_type<T_log_location, T_precision>::type
neg_binomial_2_log_lpmf(const T_n& n, const T_log_location& eta,
                        const T_precision& phi) {
  return neg_binomial_2_log_lpmf<false>(n, eta, phi);
}
}  // namespace math
}  // namespace stan
#endif
