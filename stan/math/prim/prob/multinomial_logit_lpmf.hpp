#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
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
template <bool propto, typename T_beta, typename T_prob = scalar_type_t<T_beta>,
          require_eigen_col_vector_t<T_beta>* = nullptr>
return_type_t<T_prob> multinomial_logit_lpmf(const std::vector<int>& ns,
                                             const T_beta& beta) {
  static constexpr const char* function = "multinomial_logit_lpmf";
  check_size_match(function, "Size of number of trials variable", ns.size(),
                   "rows of log-probabilities parameter", beta.rows());
  check_nonnegative(function, "Number of trials variable", ns);
  const auto& beta_ref = to_ref(beta);
  check_finite(function, "log-probabilities parameter", beta_ref);

  return_type_t<T_prob> lp(0.0);

  decltype(auto) ns_map = as_array_or_scalar(ns);

  if (include_summand<propto>::value) {
    lp += lgamma(1 + ns_map.sum()) - lgamma(1 + ns_map).sum();
  }

  if (include_summand<propto, T_prob>::value) {
    T_prob alpha = log_sum_exp(beta_ref);
    for (unsigned int i = 0; i < ns.size(); ++i) {
      if (ns[i] != 0) {
        lp += ns[i] * (beta_ref.coeff(i) - alpha);
      }
    }
  }

  return lp;
}

template <typename T_beta, require_eigen_col_vector_t<T_beta>* = nullptr>
return_type_t<T_beta> multinomial_logit_lpmf(const std::vector<int>& ns,
                                             const T_beta& beta) {
  return multinomial_logit_lpmf<false>(ns, beta);
}

}  // namespace math
}  // namespace stan
#endif
