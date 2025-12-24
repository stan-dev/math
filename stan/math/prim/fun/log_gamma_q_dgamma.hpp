#ifndef STAN_MATH_PRIM_FUN_LOG_GAMMA_Q_DGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LOG_GAMMA_Q_DGAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Result structure containing log(Q(a,z)) and its gradient with respect to a.
 *
 * @tparam T return type
 */
template <typename T>
struct log_gamma_q_result {
  T log_q;      ///< log(Q(a,z)) where Q is upper regularized incomplete gamma
  T dlog_q_da;  ///< d/da log(Q(a,z))
};

/**
 * Compute log(Q(a,z)) and its gradient with respect to a using continued
 * fraction expansion, where Q(a,z) = Gamma(a,z) / Gamma(a) is the regularized
 * upper incomplete gamma function.
 *
 * This uses a continued fraction representation for numerical stability when
 * computing the upper incomplete gamma function in log space, along with
 * analytical gradient computation.
 *
 * @tparam T_a type of the shape parameter
 * @tparam T_z type of the value parameter
 * @param a shape parameter (must be positive)
 * @param z value parameter (must be non-negative)
 * @param max_steps maximum iterations for continued fraction
 * @param precision convergence threshold
 * @return structure containing log(Q(a,z)) and d/da log(Q(a,z))
 */
template <typename T_a, typename T_z>
inline log_gamma_q_result<return_type_t<T_a, T_z>> log_gamma_q_dgamma(
    const T_a& a, const T_z& z, int max_steps = 250,
    double precision = 1e-16) {
  using std::exp;
  using std::fabs;
  using std::log;
  using T_return = return_type_t<T_a, T_z>;

  const double a_dbl = value_of(a);
  const double z_dbl = value_of(z);

  log_gamma_q_result<T_return> result;

  // For z > a + 1, use continued fraction for better numerical stability
  if (z_dbl > a_dbl + 1.0) {
    // Continued fraction for Q(a,z) in log space
    // log(Q(a,z)) = log_prefactor - log(continued_fraction)
    const double log_prefactor = a_dbl * log(z_dbl) - z_dbl - lgamma(a_dbl);

    double b = z_dbl + 1.0 - a_dbl;
    double C = (fabs(b) >= EPSILON) ? b : EPSILON;
    double D = 0.0;
    double f = C;

    for (int i = 1; i <= max_steps; ++i) {
      const double an = -i * (i - a_dbl);
      b += 2.0;

      D = b + an * D;
      if (fabs(D) < EPSILON) {
        D = EPSILON;
      }
      C = b + an / C;
      if (fabs(C) < EPSILON) {
        C = EPSILON;
      }

      D = 1.0 / D;
      const double delta = C * D;
      f *= delta;

      const double delta_m1 = fabs(delta - 1.0);
      if (delta_m1 < precision) {
        break;
      }
    }

    result.log_q = log_prefactor - log(f);

    // For gradient, use: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const double Q_val = exp(result.log_q);
    const double dQ_da = grad_reg_inc_gamma(a_dbl, z_dbl, tgamma(a_dbl), digamma(a_dbl));
    result.dlog_q_da = dQ_da / Q_val;

  } else {
    // For z <= a + 1, use log1m(P(a,z)) for better numerical accuracy
    const double P_val = gamma_p(a_dbl, z_dbl);
    result.log_q = log1m(P_val);

    // Gradient: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const double Q_val = exp(result.log_q);
    if (Q_val > 0) {
      const double dQ_da = grad_reg_inc_gamma(a_dbl, z_dbl, tgamma(a_dbl), digamma(a_dbl));
      result.dlog_q_da = dQ_da / Q_val;
    } else {
      // Fallback if Q rounds to zero - use asymptotic approximation
      result.dlog_q_da = log(z_dbl) - digamma(a_dbl);
    }
  }

  return result;
}

}  // namespace math
}  // namespace stan

#endif
