#ifndef STAN_MATH_PRIM_PROB_BINOMIAL_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BINOMIAL_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/log_inv_logit.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

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
template <bool propto, typename T_n, typename T_N, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_N, T_prob>* = nullptr>
return_type_t<T_prob> binomial_logit_lpmf(const T_n& n, const T_N& N,
                                          const T_prob& alpha) {
  using T_partials_return = partials_return_t<T_n, T_N, T_prob>;
  using T_n_ref = ref_type_if_not_constant_t<T_n>;
  using T_N_ref = ref_type_if_not_constant_t<T_N>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_prob>;
  static const char* function = "binomial_logit_lpmf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", alpha);

  T_n_ref n_ref = n;
  T_N_ref N_ref = N;
  T_alpha_ref alpha_ref = alpha;

  decltype(auto) n_val = to_ref(as_value_column_array_or_scalar(n_ref));
  decltype(auto) N_val = to_ref(as_value_column_array_or_scalar(N_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_bounded(function, "Successes variable", value_of(n_val), 0, N_val);
  check_nonnegative(function, "Population size parameter", N_val);
  check_finite(function, "Probability parameter", alpha_val);

  if (size_zero(n, N, alpha)) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }
  const auto& log_inv_logit_alpha
      = to_ref_if<!is_constant_all<T_prob>::value>(log_inv_logit(alpha_val));
  const auto& log1m_inv_logit_alpha
      = to_ref_if<!is_constant_all<T_prob>::value>(log1m_inv_logit(alpha_val));

  size_t maximum_size = max_size(n, N, alpha);
  T_partials_return logp = sum(n_val * log_inv_logit_alpha
                               + (N_val - n_val) * log1m_inv_logit_alpha);
  if (include_summand<propto, T_n, T_N>::value) {
    logp += sum(binomial_coefficient_log(N_val, n_val)) * maximum_size
            / max_size(n, N);
  }

  auto ops_partials = make_partials_propagator(alpha_ref);
  if (!is_constant_all<T_prob>::value) {
    edge<0>(ops_partials).partials_ = n_val - N_val * exp(log_inv_logit_alpha);
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
