#ifndef STAN_MATH_PRIM_PROB_BINOMIAL_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BINOMIAL_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Binomial log PMF in logit parametrization. Binomial(n|n, inv_logit(alpha))
 *
 * If given vectors of matching lengths, returns
 * the log sum of probabilities.
 *
 * @param n successes variable
 * @param N population size parameter
 * @param alpha logit transformed probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if N is negative or probability parameter is invalid
 * @throw std::invalid_argument if vector sizes do not match
 */
template <bool propto, typename T_n, typename T_N, typename T_prob>
return_type_t<T_prob> binomial_logit_lpmf(const T_n& n, const T_N& N,
                                          const T_prob& alpha) {
  using T_partials_return = partials_return_t<T_n, T_N, T_prob>;

  static const char* function = "binomial_logit_lpmf";

  if (size_zero(n, N, alpha)) {
    return 0.0;
  }

  T_partials_return logp = 0;
  check_bounded(function, "Successes variable", n, 0, N);
  check_nonnegative(function, "Population size parameter", N);
  check_finite(function, "Probability parameter", alpha);
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", alpha);

  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  using std::log;

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_N> N_vec(N);
  scalar_seq_view<T_prob> alpha_vec(alpha);
  size_t size_alpha = stan::math::size(alpha);
  size_t max_size_seq_view = max_size(n, N, alpha);

  operands_and_partials<T_prob> ops_partials(alpha);

  if (include_summand<propto>::value) {
    for (size_t i = 0; i < max_size_seq_view; ++i) {
      logp += binomial_coefficient_log(N_vec[i], n_vec[i]);
    }
  }

  VectorBuilder<true, T_partials_return, T_prob> inv_logit_alpha(size_alpha);
  VectorBuilder<true, T_partials_return, T_prob> inv_logit_neg_alpha(
      size_alpha);
  VectorBuilder<true, T_partials_return, T_prob> log_inv_logit_alpha(
      size_alpha);
  VectorBuilder<true, T_partials_return, T_prob> log_inv_logit_neg_alpha(
      size_alpha);

  for (size_t i = 0; i < size_alpha; ++i) {
    const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
    inv_logit_alpha[i] = inv_logit(alpha_dbl);
    inv_logit_neg_alpha[i] = inv_logit(-alpha_dbl);
    log_inv_logit_alpha[i] = log(inv_logit_alpha[i]);
    log_inv_logit_neg_alpha[i] = log(inv_logit_neg_alpha[i]);
  }

  for (size_t i = 0; i < max_size_seq_view; ++i) {
    logp += n_vec[i] * log_inv_logit_alpha[i]
            + (N_vec[i] - n_vec[i]) * log_inv_logit_neg_alpha[i];
  }

  if (!is_constant_all<T_prob>::value) {
    if (size_alpha == 1) {
      T_partials_return sum_n = 0;
      T_partials_return sum_N = 0;
      for (size_t i = 0; i < max_size_seq_view; ++i) {
        sum_n += n_vec[i];
        sum_N += N_vec[i];
      }
      ops_partials.edge1_.partials_[0]
          += sum_n * inv_logit_neg_alpha[0]
             - (sum_N - sum_n) * inv_logit_alpha[0];
    } else {
      for (size_t i = 0; i < max_size_seq_view; ++i) {
        ops_partials.edge1_.partials_[i]
            += n_vec[i] * inv_logit_neg_alpha[i]
               - (N_vec[i] - n_vec[i]) * inv_logit_alpha[i];
      }
    }
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_N, typename T_prob>
inline return_type_t<T_prob> binomial_logit_lpmf(const T_n& n, const T_N& N,
                                                 const T_prob& alpha) {
  return binomial_logit_lpmf<false>(n, N, alpha);
}

}  // namespace math
}  // namespace stan
#endif
