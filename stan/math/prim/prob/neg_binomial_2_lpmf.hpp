#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
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
  size_t max_size_seq_view = max_size(n, mu, phi);

  operands_and_partials<T_location, T_precision> ops_partials(mu, phi);

  size_t len_ep = max_size(mu, phi);
  size_t len_np = max_size(n, phi);

  VectorBuilder<true, T_partials_return, T_location> mu__(size(mu));

  for (size_t i = 0, max_size_seq_view = size(mu); i < max_size_seq_view; ++i) {
    mu__[i] = value_of(mu_vec[i]);
  }

  VectorBuilder<true, T_partials_return, T_precision> phi__(size(phi));
  for (size_t i = 0, max_size_seq_view = size(phi); i < max_size_seq_view;
       ++i) {
    phi__[i] = value_of(phi_vec[i]);
  }

  VectorBuilder<true, T_partials_return, T_precision> log_phi(size(phi));
  for (size_t i = 0, max_size_seq_view = size(phi); i < max_size_seq_view;
       ++i) {
    log_phi[i] = log(phi__[i]);
  }

  VectorBuilder<true, T_partials_return, T_location, T_precision>
      log_mu_plus_phi(len_ep);
  for (size_t i = 0; i < len_ep; ++i) {
    log_mu_plus_phi[i] = log(mu__[i] + phi__[i]);
  }

  VectorBuilder<true, T_partials_return, T_n, T_precision> n_plus_phi(len_np);
  for (size_t i = 0; i < len_np; ++i) {
    n_plus_phi[i] = n_vec[i] + phi__[i];
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // if phi is large we probably overflow, defer to Poisson:
    if (phi__[i] > 1e5) {
      // TODO(martinmodrak) This is wrong (doesn't pass propto information),
      // and inaccurate for n = 0, but shouldn't break most models.
      // Also the 1e5 cutoff is too small.
      // Will be adressed better in PR #1497
      logp += poisson_lpmf(n_vec[i], mu__[i]);
    } else {
      if (include_summand<propto>::value) {
        logp -= lgamma(n_vec[i] + 1.0);
      }
      if (include_summand<propto, T_precision>::value) {
        logp += multiply_log(phi__[i], phi__[i]) - lgamma(phi__[i]);
      }
      if (include_summand<propto, T_location>::value) {
        logp += multiply_log(n_vec[i], mu__[i]);
      }
      if (include_summand<propto, T_precision>::value) {
        logp += lgamma(n_plus_phi[i]);
      }
      logp -= (n_plus_phi[i]) * log_mu_plus_phi[i];
    }

    if (!is_constant_all<T_location>::value) {
      ops_partials.edge1_.partials_[i]
          += n_vec[i] / mu__[i] - (n_vec[i] + phi__[i]) / (mu__[i] + phi__[i]);
    }
    if (!is_constant_all<T_precision>::value) {
      ops_partials.edge2_.partials_[i]
          += 1.0 - n_plus_phi[i] / (mu__[i] + phi__[i]) + log_phi[i]
             - log_mu_plus_phi[i] - digamma(phi__[i]) + digamma(n_plus_phi[i]);
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
