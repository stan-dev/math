#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_CDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/inc_beta.hpp>
#include <stan/math/prim/scal/fun/inc_beta_dda.hpp>
#include <stan/math/prim/scal/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_location, typename T_precision>
inline auto neg_binomial_2_cdf(T_n&& n, T_location&& mu,
                               T_precision&& phi) {
  using T_partials = partials_return_t<T_n, T_location, T_precision>;
  T_partials P(1.0);


  static const char* function = "neg_binomial_2_cdf";
  check_positive_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Precision parameter", phi);
  check_not_nan(function, "Random variable", n);
  check_consistent_sizes(function, "Random variable", n, "Location parameter",
                         mu, "Precision Parameter", phi);

  const scalar_seq_view<T_n> n_vec(n);
  const scalar_seq_view<T_location> mu_vec(mu);
  const scalar_seq_view<T_precision> phi_vec(phi);
  const size_t size = max_size(n, mu, phi);
  operands_and_partials<T_location, T_precision> ops_partials(mu, phi);

  if (size_zero(n, mu, phi)) {
    return ops_partials.build(P);
  }

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::length(n); i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(T_partials(0.0));
    }
  }

  VectorBuilder<!is_constant_all<T_precision>::value, T_partials, T_precision>
      digamma_phi_vec(stan::length(phi));

  VectorBuilder<!is_constant_all<T_precision>::value, T_partials, T_precision>
      digamma_sum_vec(stan::length(phi));

  if (!is_constant_all<T_precision>::value) {
    for (size_t i = 0; i < stan::length(phi); i++) {
      const T_partials n_dbl = value_of(n_vec[i]);
      const T_partials phi_dbl = value_of(phi_vec[i]);

      digamma_phi_vec[i] = digamma(phi_dbl);
      digamma_sum_vec[i] = digamma(n_dbl + phi_dbl + 1);
    }
  }

  for (size_t i = 0; i < size; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) == std::numeric_limits<int>::max()) {
      return ops_partials.build(T_partials(1.0));
    }

    const T_partials n_dbl = value_of(n_vec[i]);
    const T_partials mu_dbl = value_of(mu_vec[i]);
    const T_partials phi_dbl = value_of(phi_vec[i]);

    const T_partials p_dbl = phi_dbl / (mu_dbl + phi_dbl);
    const T_partials d_dbl = 1.0 / ((mu_dbl + phi_dbl) * (mu_dbl + phi_dbl));

    const T_partials P_i = inc_beta(phi_dbl, n_dbl + 1.0, p_dbl);

    P *= P_i;

    if (!is_constant_all<T_location>::value) {
      ops_partials.edge1_.partials_[i]
          += -inc_beta_ddz(phi_dbl, n_dbl + 1.0, p_dbl) * phi_dbl * d_dbl / P_i;
    }

    if (!is_constant_all<T_precision>::value) {
      ops_partials.edge2_.partials_[i]
          += inc_beta_dda(phi_dbl, n_dbl + 1, p_dbl, digamma_phi_vec[i],
                          digamma_sum_vec[i])
                 / P_i
             + inc_beta_ddz(phi_dbl, n_dbl + 1.0, p_dbl) * mu_dbl * d_dbl / P_i;
    }
  }

  if (!is_constant_all<T_location>::value) {
    for (size_t i = 0; i < stan::length(mu); ++i) {
      ops_partials.edge1_.partials_[i] *= P;
    }
  }

  if (!is_constant_all<T_precision>::value) {
    for (size_t i = 0; i < stan::length(phi); ++i) {
      ops_partials.edge2_.partials_[i] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
