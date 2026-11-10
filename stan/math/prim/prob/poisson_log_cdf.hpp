#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_CDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log_gamma_q_dgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The cumulative distribution function of the Poisson distribution with
 * log rate alpha, Q(n + 1, exp(alpha)) with the regularized upper
 * incomplete gamma function Q. It is computed from log Q, so its gradient
 * stays finite where the probability underflows.
 *
 * @tparam T_n type of n
 * @tparam T_log_rate type of alpha
 * @param n random variable
 * @param alpha log rate parameter
 * @return probability, or product of probabilities for containers
 * @throw std::domain_error if alpha is NaN
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_log_rate>
inline return_type_t<T_log_rate> poisson_log_cdf(const T_n& n,
                                                 const T_log_rate& alpha) {
  using T_partials_return = partials_return_t<T_n, T_log_rate>;
  using T_n_ref = ref_type_if_not_constant_t<T_n>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_log_rate>;
  static constexpr const char* function = "poisson_log_cdf";
  check_consistent_sizes(function, "Random variable", n, "Log rate parameter",
                         alpha);
  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  check_not_nan(function, "Log rate parameter", alpha_ref);

  if (size_zero(n, alpha)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(alpha_ref);
  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  const size_t N = max_size(n, alpha);

  // An element with probability 0 (n < 0, or a rate exp(alpha) that is
  // infinite) makes the result 0; the gradients are treated as zero
  for (size_t i = 0; i < N; i++) {
    if (n_vec.val(i) < 0 || std::exp(value_of_rec(alpha_vec.val(i))) == INFTY) {
      return ops_partials.build(0.0);
    }
  }

  T_partials_return log_P(0.0);
  for (size_t i = 0; i < N; i++) {
    const T_partials_return lambda = exp(T_partials_return(alpha_vec.val(i)));
    // An element with probability 1 (n = INT_MAX, which stands for an
    // infinite count, or a rate of 0) contributes nothing
    if (n_vec.val(i) == std::numeric_limits<int>::max()
        || value_of_rec(lambda) == 0) {
      continue;
    }
    const auto tail
        = internal::log_gamma_tail<true, false, is_autodiff_v<T_log_rate>>(
            n_vec.val(i) + 1.0, lambda);
    log_P += tail.log_prob;
    if constexpr (is_autodiff_v<T_log_rate>) {
      partials<0>(ops_partials)[i] += tail.d_z * lambda;
    }
  }
  const T_partials_return P = exp(log_P);
  if constexpr (is_autodiff_v<T_log_rate>) {
    for (size_t i = 0; i < stan::math::size(alpha); ++i) {
      partials<0>(ops_partials)[i] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
