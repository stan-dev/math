#ifndef STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/beta_binomial_log_cdfs.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the Beta-Binomial distribution with given
 * population size, prior success, and prior failure parameters. Given
 * containers of matching sizes, returns the log sum of probabilities.
 *
 * The lower or the upper tail is summed directly, and the other one is its
 * complement only where that loses no digits; see
 * internal::beta_binomial_log_cdfs().
 *
 * @tparam T_n type of success parameter
 * @tparam T_N type of population size parameter
 * @tparam T_size1 type of prior success parameter
 * @tparam T_size2 type of prior failure parameter
 *
 * @param n success parameter
 * @param N population size parameter
 * @param alpha prior success parameter
 * @param beta prior failure parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if N, alpha, or beta fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_N, typename T_size1, typename T_size2>
inline return_type_t<T_size1, T_size2> beta_binomial_lccdf(
    const T_n& n, const T_N& N, const T_size1& alpha, const T_size2& beta) {
  using T_partials_return = partials_return_t<T_n, T_N, T_size1, T_size2>;
  using T_N_ref = ref_type_t<T_N>;
  using T_alpha_ref = ref_type_t<T_size1>;
  using T_beta_ref = ref_type_t<T_size2>;
  static constexpr const char* function = "beta_binomial_lccdf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "First prior sample size parameter", alpha,
                         "Second prior sample size parameter", beta);
  if (size_zero(n, N, alpha, beta)) {
    return 0;
  }

  T_N_ref N_ref = N;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "Population size parameter", N_ref);
  check_positive_finite(function, "First prior sample size parameter",
                        alpha_ref);
  check_positive_finite(function, "Second prior sample size parameter",
                        beta_ref);

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(alpha_ref, beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_N_ref> N_vec(N_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t max_size_seq_view = max_size(n, N, alpha, beta);

  constexpr bool any_autodiff = is_any_autodiff_v<T_size1, T_size2>;
  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_vec.val(i) < 0) {
      continue;  // a term of log(1) = 0
    }
    if (n_vec.val(i) >= N_vec.val(i)) {
      return ops_partials.build(negative_infinity());
    }
    const auto res = internal::beta_binomial_log_cdfs<any_autodiff>(
        n_vec.val(i), N_vec.val(i), T_partials_return(alpha_vec.val(i)),
        T_partials_return(beta_vec.val(i)));
    P += res.lccdf;
    if constexpr (is_autodiff_v<T_size1>) {
      partials<0>(ops_partials)[i] += res.dlccdf[0];
    }
    if constexpr (is_autodiff_v<T_size2>) {
      partials<1>(ops_partials)[i] += res.dlccdf[1];
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
