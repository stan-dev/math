#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/prob/neg_binomial_ccdf_log.hpp>

namespace stan {
namespace math {

// Temporary neg_binomial_2_ccdf implementation that
// transforms the input parameters and calls neg_binomial_ccdf
template <typename T_n, typename T_location, typename T_precision>
return_type_t<T_location, T_precision> neg_binomial_2_lccdf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  using T_mu_ref = ref_type_t<T_location>;
  using T_phi_ref = ref_type_t<T_precision>;
  static constexpr const char* function = "neg_binomial_2_lccdf";
  check_consistent_sizes(function, "Random variable", n, "Location parameter",
                         mu, "Precision Parameter", phi);
  T_mu_ref mu_ref = mu;
  T_phi_ref phi_ref = phi;
  check_positive_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Precision parameter", phi_ref);
  check_not_nan(function, "Random variable", n);

  if (size_zero(n, mu, phi)) {
    return 0;
  }

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t size_beta = max_size(mu, phi);

  VectorBuilder<true, return_type_t<T_location, T_precision>, T_location,
                T_precision>
      beta_vec(size_beta);
  for (size_t i = 0; i < size_beta; ++i) {
    beta_vec[i] = phi_vec[i] / mu_vec[i];
  }

  return neg_binomial_lccdf(n, phi_ref, beta_vec.data());
}

}  // namespace math
}  // namespace stan
#endif
