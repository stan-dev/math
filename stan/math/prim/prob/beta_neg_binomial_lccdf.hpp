#ifndef STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/beta_neg_binomial_log_cdfs.hpp>
#include <limits>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the Beta-Negative Binomial distribution with given
 * number of successes, prior success, and prior failure parameters.
 * Given containers of matching sizes, returns the log sum of probabilities.
 *
 * The lower or the upper tail is summed directly, and the other one is its
 * complement; see internal::beta_neg_binomial_log_cdfs().
 *
 * @tparam T_n type of failure parameter
 * @tparam T_r type of number of successes parameter
 * @tparam T_alpha type of prior success parameter
 * @tparam T_beta type of prior failure parameter
 *
 * @param n failure parameter
 * @param r Number of successes parameter
 * @param alpha prior success parameter
 * @param beta prior failure parameter
 * @param precision not used; the sums stop when their remaining terms are
 *   below 1e-17 of the partial sums
 * @param max_steps largest number of terms of a sum, default \f$10^{8}\f$
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if r, alpha, or beta fails to be positive, or if
 *   a sum needs more than max_steps terms
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_r, typename T_alpha, typename T_beta>
inline return_type_t<T_r, T_alpha, T_beta> beta_neg_binomial_lccdf(
    const T_n& n, const T_r& r, const T_alpha& alpha, const T_beta& beta,
    const double precision = 1e-8, const int max_steps = 1e8) {
  static constexpr const char* function = "beta_neg_binomial_lccdf";
  check_consistent_sizes(
      function, "Failures variable", n, "Number of successes parameter", r,
      "Prior success parameter", alpha, "Prior failure parameter", beta);
  if (size_zero(n, r, alpha, beta)) {
    return 0;
  }

  using T_r_ref = ref_type_t<T_r>;
  T_r_ref r_ref = r;
  using T_alpha_ref = ref_type_t<T_alpha>;
  T_alpha_ref alpha_ref = alpha;
  using T_beta_ref = ref_type_t<T_beta>;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Number of successes parameter", r_ref);
  check_positive_finite(function, "Prior success parameter", alpha_ref);
  check_positive_finite(function, "Prior failure parameter", beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_r_ref> r_vec(r_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  int size_n = stan::math::size(n);
  size_t max_size_seq_view = max_size(n, r, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (int i = 0; i < size_n; i++) {
    if (n_vec.val(i) < 0) {
      return 0.0;
    }
  }

  using T_partials_return = partials_return_t<T_n, T_r, T_alpha, T_beta>;
  constexpr bool any_autodiff = is_any_autodiff_v<T_r, T_alpha, T_beta>;
  T_partials_return log_ccdf(0.0);
  auto ops_partials = make_partials_propagator(r_ref, alpha_ref, beta_ref);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit return for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return ops_partials.build(negative_infinity());
    }
    const auto res = internal::beta_neg_binomial_log_cdfs<any_autodiff>(
        n_vec.val(i), T_partials_return(r_vec.val(i)),
        T_partials_return(alpha_vec.val(i)), T_partials_return(beta_vec.val(i)),
        max_steps, function);
    log_ccdf += res.lccdf;
    if constexpr (is_autodiff_v<T_r>) {
      partials<0>(ops_partials)[i] += res.dlccdf[0];
    }
    if constexpr (is_autodiff_v<T_alpha>) {
      partials<1>(ops_partials)[i] += res.dlccdf[1];
    }
    if constexpr (is_autodiff_v<T_beta>) {
      partials<2>(ops_partials)[i] += res.dlccdf[2];
    }
  }

  return ops_partials.build(log_ccdf);
}

}  // namespace math
}  // namespace stan
#endif
