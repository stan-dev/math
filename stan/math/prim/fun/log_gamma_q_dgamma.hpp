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

namespace internal {

/**
 * Compute log(Q(a,z)) using continued fraction expansion for upper incomplete
 * gamma function.
 *
 * @tparam T_a Type of shape parameter a (double or fvar types)
 * @tparam T_z Type of value parameter z (double or fvar types)
 * @param a Shape parameter
 * @param z Value at which to evaluate
 * @param max_steps Maximum number of continued fraction iterations
 * @param precision Convergence threshold
 * @return log(Q(a,z)) with same type as T_a and T_z
 */
template <typename T_a, typename T_z>
inline auto log_q_gamma_cf(const T_a& a, const T_z& z, int max_steps = 250,
                           double precision = 1e-16) {
  using stan::math::lgamma;
  using stan::math::log;
  using stan::math::value_of;
  using std::fabs;
  using T_return = return_type_t<T_a, T_z>;

  const T_return a_ret = a;
  const T_return z_ret = z;
  const auto log_prefactor = a_ret * log(z_ret) - z_ret - lgamma(a_ret);

  auto b = z_ret + 1.0 - a_ret;
  auto C = (fabs(value_of(b)) >= EPSILON) ? b : T_return(EPSILON);
  auto D = T_return(0.0);
  auto f = C;

  for (int i = 1; i <= max_steps; ++i) {
    auto an = -i * (i - a_ret);
    b += 2.0;

    D = b + an * D;
    if (fabs(value_of(D)) < EPSILON) {
      D = T_return(EPSILON);
    }
    C = b + an / C;
    if (fabs(value_of(C)) < EPSILON) {
      C = T_return(EPSILON);
    }

    D = 1.0 / D;
    auto delta = C * D;
    f *= delta;

    const double delta_m1 = value_of(fabs(value_of(delta) - 1.0));
    if (delta_m1 < precision) {
      break;
    }
  }

  return log_prefactor - log(f);
}

}  // namespace internal

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
    const T_a& a, const T_z& z, int max_steps = 250, double precision = 1e-16) {
  using std::exp;
  using std::log;
  using T_return = return_type_t<T_a, T_z>;

  const double a_dbl = value_of(a);
  const double z_dbl = value_of(z);

  log_gamma_q_result<T_return> result;

  // For z > a + 1, use continued fraction for better numerical stability
  if (z_dbl > a_dbl + 1.0) {
    result.log_q = internal::log_q_gamma_cf(a_dbl, z_dbl, max_steps, precision);

    // For gradient, use: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const double Q_val = exp(result.log_q);
    const double dQ_da
        = grad_reg_inc_gamma(a_dbl, z_dbl, tgamma(a_dbl), digamma(a_dbl));
    result.dlog_q_da = dQ_da / Q_val;

  } else {
    // For z <= a + 1, use log1m(P(a,z)) for better numerical accuracy
    const double P_val = gamma_p(a_dbl, z_dbl);
    result.log_q = log1m(P_val);

    // Gradient: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const double Q_val = exp(result.log_q);
    if (Q_val > 0) {
      const double dQ_da
          = grad_reg_inc_gamma(a_dbl, z_dbl, tgamma(a_dbl), digamma(a_dbl));
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
