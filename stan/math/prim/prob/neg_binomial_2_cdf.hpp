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
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_location, typename T_precision>
return_type_t<T_location, T_precision> neg_binomial_2_cdf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  using T_partials_return = partials_return_t<T_n, T_location, T_precision>;
  using T_n_ref = ref_type_t<T_n>;
  using T_mu_ref = ref_type_t<T_location>;
  using T_phi_ref = ref_type_t<T_precision>;
  static constexpr const char* function = "neg_binomial_2_cdf";
  check_consistent_sizes(function, "Random variable", n, "Location parameter",
                         mu, "Precision Parameter", phi);

  T_n_ref n_ref = n;
  T_mu_ref mu_ref = mu;
  T_phi_ref phi_ref = phi;

  check_positive_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Precision parameter", phi_ref);
  check_not_nan(function, "Random variable", n_ref);

  if (size_zero(n, mu, phi)) {
    return 1.0;
  }

  T_partials_return P(1.0);
  auto ops_partials = make_partials_propagator(mu_ref, phi_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t size_phi = stan::math::size(phi);
  size_t size_n_phi = max_size(n, phi);
  size_t max_size_seq_view = max_size(n, mu, phi);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 0) {
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
      digamma_phi_vec[i] = digamma(phi_vec.val(i));
    }
    for (size_t i = 0; i < size_n_phi; i++) {
      const T_partials_return n_dbl = n_vec.val(i);
      const T_partials_return phi_dbl = phi_vec.val(i);
      digamma_sum_vec[i] = digamma(n_dbl + phi_dbl + 1);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(1.0);
    }

    const T_partials_return n_dbl_p1 = n_vec.val(i) + 1;
    const T_partials_return mu_dbl = mu_vec.val(i);
    const T_partials_return phi_dbl = phi_vec.val(i);
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
      partials<0>(ops_partials)[i] -= inc_beta_ddz_i * phi_dbl;
    }

    if (!is_constant_all<T_precision>::value) {
      partials<1>(ops_partials)[i]
          += inc_beta_dda(phi_dbl, n_dbl_p1, p_dbl, digamma_phi_vec[i],
                          digamma_sum_vec[i])
                 / P_i
             + inc_beta_ddz_i * mu_dbl;
    }
  }

  if (!is_constant_all<T_location>::value) {
    for (size_t i = 0; i < stan::math::size(mu); ++i) {
      partials<0>(ops_partials)[i] *= P;
    }
  }

  if (!is_constant_all<T_precision>::value) {
    for (size_t i = 0; i < size_phi; ++i) {
      partials<1>(ops_partials)[i] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
