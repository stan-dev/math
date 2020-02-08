#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_CDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_location, typename T_precision>
return_type_t<T_location, T_precision> neg_binomial_2_cdf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  static const char* function = "neg_binomial_2_cdf";
  using T_partials_return = partials_return_t<T_n, T_location, T_precision>;

  T_partials_return P(1.0);
  if (size_zero(n, mu, phi)) {
    return P;
  }

  check_positive_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Precision parameter", phi);
  check_not_nan(function, "Random variable", n);
  check_consistent_sizes(function, "Random variable", n, "Location parameter",
                         mu, "Precision Parameter", phi);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_location> mu_vec(mu);
  scalar_seq_view<T_precision> phi_vec(phi);
  size_t size_phi = size(phi);
  size_t size_n_phi = max_size(n, phi);
  size_t max_size_seq_view = max_size(n, mu, phi);

  operands_and_partials<T_location, T_precision> ops_partials(mu, phi);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(n); i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_precision>::value, T_partials_return,
                T_precision>
      digamma_phi_vec(size_phi);
  VectorBuilder<!is_constant_all<T_precision>::value, T_partials_return, T_n,
                T_precision>
      digamma_sum_vec(size_n_phi);

  if (!is_constant_all<T_precision>::value) {
    for (size_t i = 0; i < size_phi; i++) {
      digamma_phi_vec[i] = digamma(value_of(phi_vec[i]));
    }
    for (size_t i = 0; i < size_n_phi; i++) {
      const T_partials_return n_dbl = value_of(n_vec[i]);
      const T_partials_return phi_dbl = value_of(phi_vec[i]);
      digamma_sum_vec[i] = digamma(n_dbl + phi_dbl + 1);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) == std::numeric_limits<int>::max()) {
      return ops_partials.build(1.0);
    }

    const T_partials_return n_dbl_p1 = value_of(n_vec[i]) + 1;
    const T_partials_return mu_dbl = value_of(mu_vec[i]);
    const T_partials_return phi_dbl = value_of(phi_vec[i]);
    const T_partials_return inv_mu_plus_phi = inv(mu_dbl + phi_dbl);
    const T_partials_return p_dbl = phi_dbl * inv_mu_plus_phi;
    const T_partials_return d_dbl = square(inv_mu_plus_phi);
    const T_partials_return P_i = inc_beta(phi_dbl, n_dbl_p1, p_dbl);
    const T_partials_return inc_beta_ddz_i
        = is_constant_all<T_location, T_precision>::value
              ? 0
              : inc_beta_ddz(phi_dbl, n_dbl_p1, p_dbl) * d_dbl / P_i;

    P *= P_i;

    if (!is_constant_all<T_location>::value) {
      ops_partials.edge1_.partials_[i] -= inc_beta_ddz_i * phi_dbl;
    }

    if (!is_constant_all<T_precision>::value) {
      ops_partials.edge2_.partials_[i]
          += inc_beta_dda(phi_dbl, n_dbl_p1, p_dbl, digamma_phi_vec[i],
                          digamma_sum_vec[i])
                 / P_i
             + inc_beta_ddz_i * mu_dbl;
    }
  }

  if (!is_constant_all<T_location>::value) {
    for (size_t i = 0; i < size(mu); ++i) {
      ops_partials.edge1_.partials_[i] *= P;
    }
  }

  if (!is_constant_all<T_precision>::value) {
    for (size_t i = 0; i < size_phi; ++i) {
      ops_partials.edge2_.partials_[i] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
