#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LCDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CDF of the geometric distribution. Given containers of
 * matching sizes, returns the log of the product of probabilities.
 *
 * The geometric distribution has log CDF
 *
 * \f[
 *   \log P(N \le n \mid \theta)
 *     = \log\!\bigl(1 - (1 - \theta)^{n + 1}\bigr),
 *   \quad n \in \{0, 1, 2, \dots\},
 *   \quad \theta \in (0, 1].
 * \f]
 *
 * The gradient with respect to \f$\theta\f$ is
 *
 * \f[
 *   \frac{\partial}{\partial \theta} \log P(N \le n \mid \theta)
 *     = \frac{(n + 1)\,(1 - \theta)^n}{1 - (1 - \theta)^{n + 1}}.
 * \f]
 *
 * Implemented as \f$\mathrm{log1m\_exp}((n + 1)\,\mathrm{log1m}(\theta))\f$
 * for numerical stability.
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not in (0, 1]
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
inline return_type_t<T_prob> geometric_lcdf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_theta_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lcdf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  T_theta_ref theta_ref = theta;
  const auto& n_arr = as_value_column_array_or_scalar(n);
  const auto& theta_arr = as_value_column_array_or_scalar(theta_ref);
  check_positive_finite(function, "Probability parameter", theta_arr);
  check_less_or_equal(function, "Probability parameter", theta_arr, 1.0);

  if (size_zero(n, theta)) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  // log P(N <= n) = -inf for n < 0
  if (any(n_arr < 0)) {
    return ops_partials.build(NEGATIVE_INFTY);
  }

  // For theta = 1: log_q = -inf, log1m_exp(-inf) = log(1 - 0) = 0
  //   (correct: log P = 0 when success is certain).
  const auto& log1m_theta = log1m(theta_arr);
  const auto& log_q = (n_arr + 1.0) * log1m_theta;
  const auto& log_P_i = log1m_exp(log_q);
  T_partials_return logP = sum(log_P_i);

  if constexpr (is_autodiff_v<T_prob>) {
    // For n = 0: dP_dtheta = 1 (avoid 0 * log1m(1) = NaN at theta = 1).
    // For n > 0, theta = 1: dP_dtheta = 0 and P_i = 1 -> partial = 0.
    const auto& dP_dtheta = select(n_arr == 0, T_partials_return(1.0),
                                   (n_arr + 1.0) * exp(n_arr * log1m_theta));
    const auto& P_i = -expm1(log_q);
    if constexpr (is_stan_scalar_v<T_prob>) {
      partials<0>(ops_partials) = sum(elt_divide(dP_dtheta, P_i));
    } else {
      partials<0>(ops_partials) = elt_divide(dP_dtheta, P_i);
    }
  }

  return ops_partials.build(logP);
}

}  // namespace math
}  // namespace stan
#endif
