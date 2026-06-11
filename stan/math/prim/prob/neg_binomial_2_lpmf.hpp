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
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

// NegBinomial(n|mu, phi)  [mu >= 0; phi > 0;  n >= 0]
template <bool propto, typename T_n, typename T_location, typename T_precision,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_location, T_precision>* = nullptr>
inline return_type_t<T_location, T_precision> neg_binomial_2_lpmf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  using T_partials_return = partials_return_t<T_n, T_location, T_precision>;
  using std::log;
  using T_n_ref = ref_type_t<T_n>;
  using T_mu_ref = ref_type_t<T_location>;
  using T_phi_ref = ref_type_t<T_precision>;
  static constexpr const char* function = "neg_binomial_2_lpmf";
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
  if constexpr (!include_summand<propto, T_location, T_precision>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  auto ops_partials = make_partials_propagator(mu_ref, phi_ref);

  auto n_vec = as_array_or_scalar(as_column_vector_or_scalar(n_ref));
  auto mu_vec = as_array_or_scalar(as_column_vector_or_scalar(mu_ref));
  auto phi_vec = as_array_or_scalar(as_column_vector_or_scalar(phi_ref));
  decltype(auto) mu_val = value_of(mu_vec);
  decltype(auto) phi_val = value_of(phi_vec);
  auto log_phi = log(phi_val);
  auto mu_plus_phi = mu_val + phi_val;
  auto log_mu_plus_phi = log(mu_plus_phi);
  auto n_plus_phi = value_of(n_vec) + phi_val;
  constexpr bool include_precision
      = include_summand<propto, T_precision>::value;
  constexpr bool include_location = include_summand<propto, T_location>::value;
  auto logp_calc = [&]() {
    return -phi_val * (log1p(mu_val / phi_val))
           - value_of(n_vec) * log_mu_plus_phi;
  };
  if constexpr (include_precision || include_location) {
    if constexpr (include_precision && include_location) {
      logp += sum(binomial_coefficient_log(n_plus_phi - 1, n_vec)
                  + multiply_log(n_vec, mu_val) + logp_calc());
    } else if constexpr (include_precision) {
      logp
          += sum(binomial_coefficient_log(n_plus_phi - 1, n_vec) + logp_calc());
    } else if constexpr (include_location) {
      logp += sum(multiply_log(n_vec, mu_val) + logp_calc());
    }
  }
  if constexpr (is_autodiff_v<T_location>) {
    partials<0>(ops_partials)
        = n_vec / mu_val - (n_vec + phi_val) / mu_plus_phi;
  }
  if constexpr (is_autodiff_v<T_precision>) {
    auto log_term = select(mu_val < phi_val, log1p(-mu_val / mu_plus_phi),
                           log_phi - log_mu_plus_phi);
    partials<1>(ops_partials) = (mu_val - value_of(n_vec)) / mu_plus_phi
                                + log_term - digamma(phi_val)
                                + digamma(n_plus_phi);
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
