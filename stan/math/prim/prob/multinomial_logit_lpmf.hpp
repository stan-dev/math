#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Multinomial log PMF in log parametrization.
 * Multinomial(ns| softmax(beta))
 *
 * @param ns Array of outcome counts
 * @param beta Vector of unnormalized log probabilities
 * @return log probability
 */
template <bool propto, typename T_prob>
return_type_t<T_prob> multinomial_logit_lpmf(
    const std::vector<int>& ns,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& beta) {
  static const char* function = "multinomial_logit_lpmf";
  check_nonnegative(function, "Number of trials variable", ns);
  check_finite(function, "log-probabilities parameter", beta);
  check_size_match(function, "Size of number of trials variable", ns.size(),
                   "rows of log-probabilities parameter", beta.rows());

  return_type_t<T_prob> lp(0.0);

  if (include_summand<propto>::value) {
    double sum = 1.0;

    for (int n : ns) {
      sum += n;
    }
    lp += lgamma(sum);
    for (int n : ns) {
      lp -= lgamma(n + 1.0);
    }
  }

  if (include_summand<propto, T_prob>::value) {
    T_prob alpha = log_sum_exp(beta);
    for (unsigned int i = 0; i < ns.size(); ++i) {
      if (ns[i] != 0)
        lp += ns[i] * (beta[i] - alpha);
    }
  }

  return lp;
}

template <typename T_prob>
return_type_t<T_prob> multinomial_logit_lpmf(
    const std::vector<int>& ns,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& beta) {
  return multinomial_logit_lpmf<false>(ns, beta);
}

}  // namespace math
}  // namespace stan
#endif
