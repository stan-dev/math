#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

// NegBinomial(n|eta, phi)  [phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_log_location,
          typename T_precision,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_log_location, T_precision>* = nullptr>
return_type_t<T_log_location, T_precision> neg_binomial_2_log_lpmf(
    const T_n& n, const T_log_location& eta, const T_precision& phi) {
  using T_partials_return = partials_return_t<T_n, T_log_location, T_precision>;
  using std::exp;
  using std::log;
  using T_n_ref = ref_type_t<T_n>;
  using T_eta_ref = ref_type_t<T_log_location>;
  using T_phi_ref = ref_type_t<T_precision>;
  static constexpr const char* function = "neg_binomial_2_log_lpmf";
  check_consistent_sizes(function, "Failures variable", n,
                         "Log location parameter", eta, "Precision parameter",
                         phi);

  T_n_ref n_ref = n;
  T_eta_ref eta_ref = eta;
  T_phi_ref phi_ref = phi;

  check_nonnegative(function, "Failures variable", n_ref);
  check_finite(function, "Log location parameter", eta_ref);
  check_positive_finite(function, "Precision parameter", phi_ref);

  if (size_zero(n, eta, phi)) {
    return 0.0;
  }
  if (!include_summand<propto, T_log_location, T_precision>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  auto ops_partials = make_partials_propagator(eta_ref, phi_ref);

  scalar_seq_view<T_n> n_vec(n_ref);
  scalar_seq_view<T_eta_ref> eta_vec(eta_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t size_eta = stan::math::size(eta);
  size_t size_phi = stan::math::size(phi);
  size_t size_eta_phi = max_size(eta, phi);
  size_t size_n_phi = max_size(n, phi);
  size_t size_all = max_size(n, eta, phi);

  VectorBuilder<true, T_partials_return, T_log_location> eta_val(size_eta);
  for (size_t i = 0; i < size_eta; ++i) {
    eta_val[i] = eta_vec.val(i);
  }

  VectorBuilder<true, T_partials_return, T_precision> phi_val(size_phi);
  VectorBuilder<true, T_partials_return, T_precision> log_phi(size_phi);
  for (size_t i = 0; i < size_phi; ++i) {
    phi_val[i] = phi_vec.val(i);
    log_phi[i] = log(phi_val[i]);
  }

  VectorBuilder<!is_constant_all<T_log_location, T_precision>::value,
                T_partials_return, T_log_location>
      exp_eta(size_eta);
  if (!is_constant_all<T_log_location, T_precision>::value) {
    for (size_t i = 0; i < size_eta; ++i) {
      exp_eta[i] = exp(eta_val[i]);
    }
  }

  VectorBuilder<!is_constant_all<T_log_location, T_precision>::value,
                T_partials_return, T_log_location, T_precision>
      exp_eta_over_exp_eta_phi(size_eta_phi);
  if (!is_constant_all<T_log_location, T_precision>::value) {
    for (size_t i = 0; i < size_eta_phi; ++i) {
      exp_eta_over_exp_eta_phi[i] = inv(phi_val[i] / exp_eta[i] + 1);
    }
  }

  VectorBuilder<true, T_partials_return, T_log_location, T_precision>
      log1p_exp_eta_m_logphi(size_eta_phi);
  for (size_t i = 0; i < size_eta_phi; ++i) {
    log1p_exp_eta_m_logphi[i] = log1p_exp(eta_val[i] - log_phi[i]);
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
    if (include_summand<propto, T_log_location>::value) {
      logp += n_vec[i] * eta_val[i];
    }
    logp += -phi_val[i] * log1p_exp_eta_m_logphi[i]
            - n_vec[i] * (log_phi[i] + log1p_exp_eta_m_logphi[i]);

    if (!is_constant_all<T_log_location>::value) {
      partials<0>(ops_partials)[i]
          += n_vec[i] - n_plus_phi[i] * exp_eta_over_exp_eta_phi[i];
    }
    if (!is_constant_all<T_precision>::value) {
      partials<1>(ops_partials)[i]
          += exp_eta_over_exp_eta_phi[i] - n_vec[i] / (exp_eta[i] + phi_val[i])
             - log1p_exp_eta_m_logphi[i]
             - (digamma(phi_val[i]) - digamma(n_plus_phi[i]));
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
