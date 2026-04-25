#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LPMF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the geometric distribution. If containers
 * of matching sizes are supplied, returns the log sum of probabilities.
 *
 * The geometric distribution counts the number of failures before
 * the first success: P(N = n | theta) = theta * (1 - theta)^n.
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not in [0, 1]
 * @throw std::domain_error if n is negative
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
inline return_type_t<T_prob> geometric_lpmf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using T_theta_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lpmf";
  check_consistent_sizes(function, "Outcome variable", n,
                         "Success probability parameter", theta);
  T_n_ref n_ref = n;
  T_theta_ref theta_ref = theta;
  check_nonnegative(function, "Outcome variable", n_ref);
  check_bounded(function, "Success probability parameter", value_of(theta_ref),
                0.0, 1.0);

  if (size_zero(n, theta)) {
    return 0.0;
  }
  if constexpr (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  // Probability-zero event: theta = 1 with n > 0
  // (when theta = 1 success is certain, so only n = 0 has positive mass).
  // A loop check here keeps both scalar and vector paths simple and
  // avoids ambiguities in Eigen boolean array combination.
  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t max_sz = max_size(n_ref, theta_ref);
  for (size_t i = 0; i < max_sz; i++) {
    if (value_of(theta_vec[i]) == 1.0 && n_vec[i] > 0) {
      return ops_partials.build(NEGATIVE_INFTY);
    }
  }

  const auto& n_arr = as_value_column_array_or_scalar(n_ref);
  const auto& theta_arr = as_value_column_array_or_scalar(theta_ref);

  // log P = log(theta) + n * log1m(theta).
  // The select on n == 0 avoids 0 * log1m(1) = 0 * (-inf) = NaN when
  // theta = 1; the n > 0 && theta = 1 case was already handled above.
  const auto& log1m_theta = log1m(theta_arr);
  const auto& failure_term
      = select(n_arr == 0, T_partials_return(0), n_arr * log1m_theta);
  T_partials_return logp = sum(log(theta_arr) + failure_term);

  if constexpr (is_autodiff_v<T_prob>) {
    // d/dtheta log P = 1/theta - n/(1 - theta) = 1/theta + n/(theta - 1).
    // For n = 0 the failure term contributes nothing; for n > 0 we have
    // theta != 1 here, so theta - 1 != 0.
    const auto& failure_grad = select(n_arr == 0, T_partials_return(0),
                                      n_arr * inv(theta_arr - 1.0));
    partials<0>(ops_partials) = inv(theta_arr) + failure_grad;
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_prob>
inline return_type_t<T_prob> geometric_lpmf(const T_n& n, const T_prob& theta) {
  return geometric_lpmf<false>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
