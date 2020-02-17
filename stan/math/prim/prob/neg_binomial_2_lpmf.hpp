#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/prob/poisson_lpmf.hpp>
#include <cmath>

namespace stan {
namespace math {

// NegBinomial(n|mu, phi)  [mu >= 0; phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_location, typename T_precision>
return_type_t<T_location, T_precision> neg_binomial_2_lpmf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  using T_partials_return = partials_return_t<T_n, T_location, T_precision>;

  static const char* function = "neg_binomial_2_lpmf";

  if (size_zero(n, mu, phi)) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  check_nonnegative(function, "Failures variable", n);
  check_positive_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Precision parameter", phi);
  check_consistent_sizes(function, "Failures variable", n, "Location parameter",
                         mu, "Precision parameter", phi);

  if (!include_summand<propto, T_location, T_precision>::value) {
    return 0.0;
  }

  using std::log;

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_location> mu_vec(mu);
  scalar_seq_view<T_precision> phi_vec(phi);
  size_t size_mu = size(mu);
  size_t size_phi = size(phi);
  size_t size_mu_phi = max_size(mu, phi);
  size_t size_n_phi = max_size(n, phi);
  size_t max_size_seq_view = max_size(n, mu, phi);

  operands_and_partials<T_location, T_precision> ops_partials(mu, phi);

  VectorBuilder<true, T_partials_return, T_location> mu_val(size_mu);
  for (size_t i = 0; i < size_mu; ++i) {
    mu_val[i] = value_of(mu_vec[i]);
  }

  VectorBuilder<true, T_partials_return, T_precision> phi_val(size_phi);
  VectorBuilder<true, T_partials_return, T_precision> log_phi(size_phi);
  for (size_t i = 0; i < size_phi; ++i) {
    phi_val[i] = value_of(phi_vec[i]);
    log_phi[i] = log(phi_val[i]);
  }

  VectorBuilder<true, T_partials_return, T_location, T_precision> mu_plus_phi(
      size_mu_phi);
  VectorBuilder<true, T_partials_return, T_location, T_precision>
      log_mu_plus_phi(size_mu_phi);
  for (size_t i = 0; i < size_mu_phi; ++i) {
    mu_plus_phi[i] = mu_val[i] + phi_val[i];
    log_mu_plus_phi[i] = log(mu_plus_phi[i]);
  }

  VectorBuilder<true, T_partials_return, T_n, T_precision> n_plus_phi(
      size_n_phi);
  for (size_t i = 0; i < size_n_phi; ++i) {
    n_plus_phi[i] = n_vec[i] + phi_val[i];
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // if phi is large we probably overflow, defer to Poisson:
    if (phi_val[i] > 1e5) {
      // TODO(martinmodrak) This is wrong (doesn't pass propto information),
      // and inaccurate for n = 0, but shouldn't break most models.
      // Also the 1e5 cutoff is too small.
      // Will be adressed better in PR #1497
      logp += poisson_lpmf(n_vec[i], mu_val[i]);
    } else {
      if (include_summand<propto>::value) {
        logp -= lgamma(n_vec[i] + 1.0);
      }
      if (include_summand<propto, T_precision>::value) {
        logp += multiply_log(phi_val[i], phi_val[i]) - lgamma(phi_val[i]);
      }
      if (include_summand<propto, T_location>::value) {
        logp += multiply_log(n_vec[i], mu_val[i]);
      }
      if (include_summand<propto, T_precision>::value) {
        logp += lgamma(n_plus_phi[i]);
      }
      logp -= n_plus_phi[i] * log_mu_plus_phi[i];
    }

    if (!is_constant_all<T_location>::value) {
      ops_partials.edge1_.partials_[i]
          += n_vec[i] / mu_val[i] - n_plus_phi[i] / mu_plus_phi[i];
    }
    if (!is_constant_all<T_precision>::value) {
      ops_partials.edge2_.partials_[i] += 1.0 - n_plus_phi[i] / mu_plus_phi[i]
                                          + log_phi[i] - log_mu_plus_phi[i]
                                          - digamma(phi_val[i])
                                          + digamma(n_plus_phi[i]);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_location, typename T_precision>
inline return_type_t<T_location, T_precision> neg_binomial_2_lpmf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  return neg_binomial_2_lpmf<false>(n, mu, phi);
}

}  // namespace math
}  // namespace stan
#endif
