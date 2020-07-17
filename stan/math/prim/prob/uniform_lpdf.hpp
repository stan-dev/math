#ifndef STAN_MATH_PRIM_PROB_UNIFORM_LPDF_HPP
#define STAN_MATH_PRIM_PROB_UNIFORM_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of a uniform density for the given
 * y, lower, and upper bound.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{U}}(\alpha, \beta) \\
 \log (p (y \, |\, \alpha, \beta)) &=& \log \left( \frac{1}{\beta-\alpha}
 \right) \\
 &=& \log (1) - \log (\beta - \alpha) \\
 &=& -\log (\beta - \alpha) \\
 & & \mathrm{ where } \; y \in [\alpha, \beta], \log(0) \; \mathrm{otherwise}
 \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_low type of lower bound
 * @tparam T_high type of upper bound
 * @param y A scalar variable.
 * @param alpha Lower bound.
 * @param beta Upper bound.
 * @throw std::invalid_argument if the lower bound is greater than
 *    or equal to the lower bound
 */
template <bool propto, typename T_y, typename T_low, typename T_high>
return_type_t<T_y, T_low, T_high> uniform_lpdf(const T_y& y, const T_low& alpha,
                                               const T_high& beta) {
  using T_partials_return = partials_return_t<T_y, T_low, T_high>;
  using std::log;
  static const char* function = "uniform_lpdf";
  check_not_nan(function, "Random variable", y);
  check_finite(function, "Lower bound parameter", alpha);
  check_finite(function, "Upper bound parameter", beta);
  check_consistent_sizes(function, "Random variable", y,
                         "Lower bound parameter", alpha,
                         "Upper bound parameter", beta);
  check_greater(function, "Upper bound parameter", beta, alpha);

  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_low, T_high>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  operands_and_partials<T_y, T_low, T_high> ops_partials(y, alpha, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_low> alpha_vec(alpha);
  scalar_seq_view<T_high> beta_vec(beta);
  size_t N = max_size(y, alpha, beta);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    if (y_dbl < value_of(alpha_vec[n]) || y_dbl > value_of(beta_vec[n])) {
      return LOG_ZERO;
    }
  }

  VectorBuilder<include_summand<propto, T_low, T_high>::value,
                T_partials_return, T_low, T_high>
      inv_beta_minus_alpha(max_size(alpha, beta));
  for (size_t i = 0; i < max_size(alpha, beta); i++) {
    if (include_summand<propto, T_low, T_high>::value) {
      inv_beta_minus_alpha[i]
          = 1.0 / (value_of(beta_vec[i]) - value_of(alpha_vec[i]));
    }
  }

  VectorBuilder<include_summand<propto, T_low, T_high>::value,
                T_partials_return, T_low, T_high>
      log_beta_minus_alpha(max_size(alpha, beta));
  for (size_t i = 0; i < max_size(alpha, beta); i++) {
    if (include_summand<propto, T_low, T_high>::value) {
      log_beta_minus_alpha[i]
          = log(value_of(beta_vec[i]) - value_of(alpha_vec[i]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    if (include_summand<propto, T_low, T_high>::value) {
      logp -= log_beta_minus_alpha[n];
    }

    if (!is_constant_all<T_low>::value) {
      ops_partials.edge2_.partials_[n] += inv_beta_minus_alpha[n];
    }
    if (!is_constant_all<T_high>::value) {
      ops_partials.edge3_.partials_[n] -= inv_beta_minus_alpha[n];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_low, typename T_high>
inline return_type_t<T_y, T_low, T_high> uniform_lpdf(const T_y& y,
                                                      const T_low& alpha,
                                                      const T_high& beta) {
  return uniform_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
