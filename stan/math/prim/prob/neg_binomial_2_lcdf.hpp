#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/prob/beta_cdf_log.hpp>
#include <cmath>
#include <limits>
#include <vector>

namespace stan {
namespace math {

template <typename T_n, typename T_location, typename T_precision>
inline return_type_t<T_location, T_precision> neg_binomial_2_lcdf(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  using std::log;
  using T_n_ref = ref_type_t<T_n>;
  using T_mu_ref = ref_type_t<T_location>;
  using T_phi_ref = ref_type_t<T_precision>;
  static constexpr const char* function = "neg_binomial_2_lcdf";
  check_consistent_sizes(function, "Random variable", n, "Location parameter",
                         mu, "Precision Parameter", phi);

  T_n_ref n_ref = n;
  T_mu_ref mu_ref = mu;
  T_phi_ref phi_ref = phi;

  check_positive_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Precision parameter", phi_ref);
  check_not_nan(function, "Random variable", n_ref);

  if (size_zero(n, mu, phi)) {
    return 0;
  }

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t size_n = stan::math::size(n);
  size_t size_phi_mu = max_size(mu, phi);

  bool any_n_max = false;
  for (size_t i = 0; i < size_n; i++) {
    if (n_vec[i] < 0) {
      return LOG_ZERO;
    }
    any_n_max = any_n_max || n_vec[i] == std::numeric_limits<int>::max();
  }

  if (any_n_max) {
    // An element with probability 1 (n = INT_MAX, which stands for an
    // infinite count) contributes nothing. Remove such elements before the
    // beta lcdf call: the beta lcdf gradient is not accurate when the
    // second shape is 2^31.
    const size_t max_size_seq_view = max_size(n, mu, phi);
    std::vector<return_type_t<T_location, T_precision>> phi_mu_k;
    std::vector<scalar_type_t<T_precision>> phi_k;
    std::vector<double> np1_k;
    for (size_t i = 0; i < max_size_seq_view; i++) {
      if (n_vec[i] == std::numeric_limits<int>::max()) {
        continue;
      }
      phi_mu_k.push_back(phi_vec[i] / (phi_vec[i] + mu_vec[i]));
      phi_k.push_back(phi_vec[i]);
      np1_k.push_back(n_vec[i] + 1.0);
    }
    return beta_cdf_log(phi_mu_k, phi_k, np1_k);
  }

  VectorBuilder<true, return_type_t<T_location, T_precision>, T_location,
                T_precision>
      phi_mu(size_phi_mu);
  for (size_t i = 0; i < size_phi_mu; i++) {
    phi_mu[i] = phi_vec[i] / (phi_vec[i] + mu_vec[i]);
  }

  VectorBuilder<true, return_type_t<T_n>, T_n> np1(size_n);
  for (size_t i = 0; i < size_n; i++) {
    np1[i] = n_vec[i] + 1.0;
  }

  return beta_cdf_log(phi_mu.data(), phi_ref, np1.data());
}

}  // namespace math
}  // namespace stan
#endif
