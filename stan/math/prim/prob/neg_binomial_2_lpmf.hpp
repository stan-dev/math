#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

// NegBinomial(n|mu, phi)  [mu >= 0; phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_location, typename T_precision,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_location, T_precision>* = nullptr>
return_type_t<T_location, T_precision> neg_binomial_2_lpmf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  using T_partials_return = partials_return_t<T_n, T_location, T_precision>;
  using std::log;
  using T_n_ref = ref_type_t<T_n>;
  using T_mu_ref = ref_type_t<T_location>;
  using T_phi_ref = ref_type_t<T_precision>;
  static const char* function = "neg_binomial_2_lpmf";
  check_consistent_sizes(function, "Failures variable", n, "Location parameter",
                         mu, "Precision parameter", phi);

  T_n_ref n_ref = n;
  T_mu_ref mu_ref = mu;
  T_phi_ref phi_ref = phi;

  check_nonnegative(function, "Failures variable", n_ref);
  check_positive_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Precision parameter", phi_ref);

  if (size_zero(n, mu, phi)) {
    return 0.0;
  }
  if (!include_summand<propto, T_location, T_precision>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  operands_and_partials<T_mu_ref, T_phi_ref> ops_partials(mu_ref, phi_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t size_mu = stan::math::size(mu);
  size_t size_phi = stan::math::size(phi);
  size_t size_mu_phi = max_size(mu, phi);
  size_t size_n_phi = max_size(n, phi);
  size_t size_all = max_size(n, mu, phi);

  VectorBuilder<true, T_partials_return, T_location> mu_val(size_mu);
  for (size_t i = 0; i < size_mu; ++i) {
    mu_val[i] = mu_vec.val(i);
  }

  VectorBuilder<true, T_partials_return, T_precision> phi_val(size_phi);
  VectorBuilder<true, T_partials_return, T_precision> log_phi(size_phi);
  for (size_t i = 0; i < size_phi; ++i) {
    phi_val[i] = phi_vec.val(i);
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

  for (size_t i = 0; i < size_all; i++) {
    if (include_summand<propto, T_precision>::value) {
      logp += binomial_coefficient_log(n_plus_phi[i] - 1, n_vec[i]);
    }
    if (include_summand<propto, T_location>::value) {
      logp += multiply_log(n_vec[i], mu_val[i]);
    }
    logp += -phi_val[i] * (log1p(mu_val[i] / phi_val[i]))
            - n_vec[i] * log_mu_plus_phi[i];

    if (!is_constant_all<T_location>::value) {
      ops_partials.edge1_.partials_[i]
          += n_vec[i] / mu_val[i] - (n_vec[i] + phi_val[i]) / mu_plus_phi[i];
    }
    if (!is_constant_all<T_precision>::value) {
      T_partials_return log_term;
      if (mu_val[i] < phi_val[i]) {
        log_term = log1p(-mu_val[i] / mu_plus_phi[i]);
      } else {
        log_term = log_phi[i] - log_mu_plus_phi[i];
      }
      ops_partials.edge2_.partials_[i]
          += (mu_val[i] - n_vec[i]) / mu_plus_phi[i] + log_term
             - digamma(phi_val[i]) + digamma(n_plus_phi[i]);
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
