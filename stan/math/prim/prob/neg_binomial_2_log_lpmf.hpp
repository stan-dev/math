#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/prob/poisson_log_lpmf.hpp>
#include <cmath>

namespace stan {
namespace math {

// NegBinomial(n|eta, phi)  [phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_log_location,
          typename T_precision>
return_type_t<T_log_location, T_precision> neg_binomial_2_log_lpmf(
    const T_n& n, const T_log_location& eta, const T_precision& phi) {
  using T_partials_return = partials_return_t<T_n, T_log_location, T_precision>;

  static const char* function = "neg_binomial_2_log_lpmf";

  if (size_zero(n, eta, phi)) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  check_nonnegative(function, "Failures variable", n);
  check_finite(function, "Log location parameter", eta);
  check_positive_finite(function, "Precision parameter", phi);
  check_consistent_sizes(function, "Failures variable", n,
                         "Log location parameter", eta, "Precision parameter",
                         phi);

  if (!include_summand<propto, T_log_location, T_precision>::value) {
    return 0.0;
  }

  using std::exp;
  using std::log;

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_log_location> eta_vec(eta);
  scalar_seq_view<T_precision> phi_vec(phi);
  size_t size_eta = size(eta);
  size_t size_phi = size(phi);
  size_t size_eta_phi = max_size(eta, phi);
  size_t size_n_phi = max_size(n, phi);
  size_t max_size_seq_view = max_size(n, eta, phi);

  operands_and_partials<T_log_location, T_precision> ops_partials(eta, phi);

  VectorBuilder<true, T_partials_return, T_log_location> eta_val(size_eta);
  for (size_t i = 0; i < size_eta; ++i) {
    eta_val[i] = value_of(eta_vec[i]);
  }

  VectorBuilder<!is_constant_all<T_log_location, T_precision>::value,
                T_partials_return, T_log_location>
      exp_eta(size_eta);
  if (!is_constant_all<T_log_location, T_precision>::value) {
    for (size_t i = 0; i < size_eta; ++i) {
      exp_eta[i] = exp(eta_val[i]);
    }
  }

  VectorBuilder<true, T_partials_return, T_precision> phi_val(size_phi);
  VectorBuilder<true, T_partials_return, T_precision> log_phi(size_phi);
  for (size_t i = 0; i < size_phi; ++i) {
    phi_val[i] = value_of(phi_vec[i]);
    log_phi[i] = log(phi_val[i]);
  }

  VectorBuilder<true, T_partials_return, T_log_location, T_precision>
      logsumexp_eta_logphi(size_eta_phi);
  for (size_t i = 0; i < size_eta_phi; ++i) {
    logsumexp_eta_logphi[i] = log_sum_exp(eta_val[i], log_phi[i]);
  }

  VectorBuilder<true, T_partials_return, T_n, T_precision> n_plus_phi(
      size_n_phi);
  for (size_t i = 0; i < size_n_phi; ++i) {
    n_plus_phi[i] = n_vec[i] + phi_val[i];
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    if (phi_val[i] > 1e5) {
      // TODO(martinmodrak) This is wrong (doesn't pass propto information),
      // and inaccurate for n = 0, but shouldn't break most models.
      // Also the 1e5 cutoff is way too low.
      // Will be adressed better once PR #1497 is merged
      logp += poisson_log_lpmf(n_vec[i], eta_val[i]);
    } else {
      if (include_summand<propto>::value) {
        logp -= lgamma(n_vec[i] + 1.0);
      }
      if (include_summand<propto, T_precision>::value) {
        logp += multiply_log(phi_val[i], phi_val[i]) - lgamma(phi_val[i]);
      }
      if (include_summand<propto, T_log_location>::value) {
        logp += n_vec[i] * eta_val[i];
      }
      if (include_summand<propto, T_precision>::value) {
        logp += lgamma(n_plus_phi[i]);
      }
      logp -= (n_plus_phi[i]) * logsumexp_eta_logphi[i];
    }

    if (!is_constant_all<T_log_location>::value) {
      ops_partials.edge1_.partials_[i]
          += n_vec[i] - n_plus_phi[i] / (phi_val[i] / exp_eta[i] + 1);
    }
    if (!is_constant_all<T_precision>::value) {
      ops_partials.edge2_.partials_[i]
          += 1.0 - n_plus_phi[i] / (exp_eta[i] + phi_val[i]) + log_phi[i]
             - logsumexp_eta_logphi[i] - digamma(phi_val[i])
             + digamma(n_plus_phi[i]);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_log_location, typename T_precision>
inline return_type_t<T_log_location, T_precision> neg_binomial_2_log_lpmf(
    const T_n& n, const T_log_location& eta, const T_precision& phi) {
  return neg_binomial_2_log_lpmf<false>(n, eta, phi);
}

}  // namespace math
}  // namespace stan
#endif
