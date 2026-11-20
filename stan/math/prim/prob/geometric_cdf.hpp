#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF of the geometric distribution. Given containers of
 * matching sizes, returns the product of probabilities.
 *
 * The geometric distribution has CDF
 *
 * \f[
 *   P(N \le n \mid \theta) = 1 - (1 - \theta)^{n + 1},
 *   \quad n \in \{0, 1, 2, \dots\},
 *   \quad \theta \in (0, 1].
 * \f]
 *
 * The gradient with respect to \f$\theta\f$ is
 *
 * \f[
 *   \frac{\partial}{\partial \theta} P(N \le n \mid \theta)
 *     = (n + 1)\,(1 - \theta)^n.
 * \f]
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if theta is not in (0, 1]
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
inline return_type_t<T_prob> geometric_cdf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_theta_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_cdf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  T_theta_ref theta_ref = theta;
  const auto& n_arr = as_value_column_array_or_scalar(n);
  const auto& theta_arr = as_value_column_array_or_scalar(theta_ref);
  check_positive_finite(function, "Probability parameter", theta_arr);
  check_less_or_equal(function, "Probability parameter", theta_arr, 1.0);

  if (size_zero(n, theta)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  // P(N <= n) = 0 for n < 0
  if (any(n_arr < 0)) {
    return ops_partials.build(0.0);
  }

  // Compute via -expm1((n + 1) * log1m(theta)) for numerical stability.
  // For theta = 1: log1m(1) = -inf, (n+1)*-inf = -inf (n >= 0),
  //   expm1(-inf) = -1, so P_i = 1 (certain success means N <= n always).
  const auto& log1m_theta = log1m(theta_arr);
  const auto& P_i = -expm1((n_arr + 1.0) * log1m_theta);
  const T_partials_return P = prod(P_i);

  if constexpr (is_autodiff_v<T_prob>) {
    // Compute (n + 1) * exp(n * log1m(theta)) for numerical stability.
    // For n = 0: (n+1)*exp(0) = 1; the select avoids 0 * log1m(1) = NaN
    //   when theta = 1.
    // For n > 0, theta = 1: (n+1) * exp(n * -inf) = (n+1) * 0 = 0
    //   (correct: derivative vanishes once CDF saturates at 1).
    const auto& dP_dtheta = select(n_arr == 0, T_partials_return(1.0),
                                   (n_arr + 1.0) * exp(n_arr * log1m_theta));
    if constexpr (is_stan_scalar_v<T_prob>) {
      partials<0>(ops_partials) = sum(P * elt_divide(dP_dtheta, P_i));
    } else {
      partials<0>(ops_partials) = P * elt_divide(dP_dtheta, P_i);
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
