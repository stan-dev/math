#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/prob/beta_cdf_log.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_n, typename T_location, typename T_precision>
inline auto neg_binomial_2_lcdf(const T_n& n, const T_location& mu,
                                const T_precision& phi) {
  using std::log;
  using T_partials = partials_return_t<T_n, T_location, T_precision>;
  using T_return = return_type_t<T_n, T_location, T_precision>;

  static const char* function = "neg_binomial_2_lcdf";
  check_positive_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Precision parameter", phi);
  check_not_nan(function, "Random variable", n);
  check_consistent_sizes(function, "Random variable", n, "Location parameter",
                         mu, "Precision Parameter", phi);

  const scalar_seq_view<T_n> n_vec(n);
  const scalar_seq_view<T_location> mu_vec(mu);
  const scalar_seq_view<T_precision> phi_vec(phi);
  const size_t size_phi_mu = max_size(mu, phi);
  if (size_zero(n, mu, phi)) {
    return T_return(0.0);
  }

  VectorBuilder<true, return_type_t<T_location, T_precision>, T_location,
                T_precision>
      phi_mu(size_phi_mu);
  for (size_t i = 0; i < size_phi_mu; i++) {
    phi_mu[i] = phi_vec[i] / (phi_vec[i] + mu_vec[i]);
  }

  const size_t size_n = length(n);
  VectorBuilder<true, return_type_t<T_n>, T_n> np1(size_n);
  for (size_t i = 0; i < size_n; i++) {
    if (n_vec[i] < 0) {
      return T_return(log(0.0));
    } else {
      np1[i] = n_vec[i] + 1.0;
    }
  }

  return beta_cdf_log(phi_mu.data(), phi, np1.data());
}

}  // namespace math
}  // namespace stan
#endif
