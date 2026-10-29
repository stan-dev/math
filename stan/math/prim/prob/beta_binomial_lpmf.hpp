#ifndef STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma_stirling_diff.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/log_rising_factorial_ratio.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the Beta-Binomial distribution with given population
 * size, prior success, and prior failure parameters. Given containers of
 * matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_n type of success parameter
 * @tparam T_N type of population size parameter
 * @tparam T_size1 type of prior success parameter
 * @tparam T_size2 type of prior failure parameter
 * @param n success parameter
 * @param N population size parameter
 * @param alpha prior success parameter
 * @param beta prior failure parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if N, alpha, or beta fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_N, typename T_size1,
          typename T_size2,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_N, T_size1, T_size2>* = nullptr>
inline return_type_t<T_size1, T_size2> beta_binomial_lpmf(const T_n& n,
                                                          const T_N& N,
                                                          const T_size1& alpha,
                                                          const T_size2& beta) {
  using T_partials_return = partials_return_t<T_size1, T_size2>;
  using T_N_ref = ref_type_t<T_N>;
  using T_alpha_ref = ref_type_t<T_size1>;
  using T_beta_ref = ref_type_t<T_size2>;
  static constexpr const char* function = "beta_binomial_lpmf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "First prior sample size parameter", alpha,
                         "Second prior sample size parameter", beta);
  if (size_zero(n, N, alpha, beta)) {
    return 0.0;
  }

  T_N_ref N_ref = N;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "Population size parameter", N_ref);
  check_positive_finite(function, "First prior sample size parameter",
                        alpha_ref);
  check_positive_finite(function, "Second prior sample size parameter",
                        beta_ref);

  if constexpr (!include_summand<propto, T_size1, T_size2>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  auto ops_partials = make_partials_propagator(alpha_ref, beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_N_ref> N_vec(N_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_n_N = max_size(n, N);
  size_t size_alpha_beta = max_size(alpha, beta);
  size_t max_size_seq_view = max_size(n, N, alpha, beta);

  for (size_t i = 0; i < max_size_seq_view; i++) {
    if (n_vec[i] < 0 || n_vec[i] > N_vec[i]) {
      return ops_partials.build(LOG_ZERO);
    }
  }

  VectorBuilder<include_summand<propto>::value, T_partials_return, T_n, T_N>
      normalizing_constant(size_n_N);
  for (size_t i = 0; i < size_n_N; i++)
    if constexpr (include_summand<propto>::value)
      normalizing_constant[i] = binomial_coefficient_log(N_vec[i], n_vec[i]);

  // The value is lbeta(n + alpha, N - n + beta) - lbeta(alpha, beta). For
  // large shapes the two lbeta values are of the order of the shapes and
  // their plain difference keeps no correct digits from shapes near 1e15;
  // internal::log_beta_ratio forms the difference
  // without that cancellation. Its shape-only part is computed once per
  // pair of shapes.
  VectorBuilder<true, T_partials_return, T_size1, T_size2> lbeta_denominator(
      size_alpha_beta);
  for (size_t i = 0; i < size_alpha_beta; i++) {
    lbeta_denominator[i] = internal::log_beta_ratio_denominator(
        T_partials_return(alpha_vec.val(i)),
        T_partials_return(beta_vec.val(i)));
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const T_partials_return alpha_dbl = alpha_vec.val(i);
    const T_partials_return beta_dbl = beta_vec.val(i);
    const T_partials_return n_i = n_vec.val(i);
    const T_partials_return N_i = N_vec.val(i);
    const T_partials_return m_i = N_i - n_i;
    if constexpr (include_summand<propto>::value) {
      logp += normalizing_constant[i];
    }
    logp += internal::log_beta_ratio(alpha_dbl, beta_dbl, n_i, m_i,
                                     lbeta_denominator[i]);

    // Each partial is a sum of two digamma differences psi(x + k) - psi(x)
    // of size k / x; digamma_diff forms them without cancellation.
    if constexpr (is_any_autodiff_v<T_size1, T_size2>) {
      const T_partials_return digamma_diff_total
          = digamma_diff(alpha_dbl + beta_dbl, N_i);
      if constexpr (is_autodiff_v<T_size1>) {
        partials<0>(ops_partials)[i]
            += digamma_diff(alpha_dbl, n_i) - digamma_diff_total;
      }
      if constexpr (is_autodiff_v<T_size2>) {
        partials<1>(ops_partials)[i]
            += digamma_diff(beta_dbl, m_i) - digamma_diff_total;
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_N, typename T_size1, typename T_size2>
inline return_type_t<T_size1, T_size2> beta_binomial_lpmf(const T_n& n,
                                                          const T_N& N,
                                                          const T_size1& alpha,
                                                          const T_size2& beta) {
  return beta_binomial_lpmf<false>(n, N, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
