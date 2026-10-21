#ifndef STAN_MATH_PRIM_FUN_LOG_GAMMA_Q_DGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LOG_GAMMA_Q_DGAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace stan {
namespace math {

namespace internal {

constexpr double LOG_Q_GAMMA_CF_PRECISION = 1.49012e-12;

/**
 * Compute log(Q(a,z)) using continued fraction expansion for upper incomplete
 * gamma function.
 *
 * @tparam T_a Type of shape parameter a (double or fvar types)
 * @tparam T_z Type of value parameter z (double or fvar types)
 * @param a Shape parameter
 * @param z Value at which to evaluate
 * @param precision Convergence threshold, default of sqrt(machine_epsilon)
 * @param max_steps Maximum number of continued fraction iterations
 * @return log(Q(a,z)) with the return type of T_a and T_z
 */
template <typename T_a, typename T_z>
inline return_type_t<T_a, T_z> log_q_gamma_cf(const T_a& a, const T_z& z,
                                              double precision
                                              = LOG_Q_GAMMA_CF_PRECISION,
                                              int max_steps = 250) {
  using T_return = return_type_t<T_a, T_z>;
  const T_return log_prefactor = a * log(z) - z - lgamma(a);

  T_return b_init = z + 1.0 - a;
  T_return C = (fabs(value_of_rec(b_init)) >= EPSILON)
                   ? b_init
                   : std::decay_t<decltype(b_init)>(EPSILON);
  T_return D = 0.0;
  T_return f = C;
  for (int i = 1; i <= max_steps; ++i) {
    T_a an = -i * (i - a);
    const T_return b = b_init + 2.0 * i;
    D = b + an * D;
    D = (fabs(value_of_rec(D)) >= EPSILON) ? D
                                           : std::decay_t<decltype(D)>(EPSILON);
    C = b + an / C;
    C = (fabs(value_of_rec(C)) >= EPSILON) ? C
                                           : std::decay_t<decltype(C)>(EPSILON);
    D = inv(D);
    const T_return delta = C * D;
    f *= delta;
    const double delta_m1 = fabs(value_of_rec(delta) - 1.0);
    if (delta_m1 < precision) {
      if constexpr (is_fvar<T_return>::value) {
        // The tangent converges more slowly than the value near z = a + 1.
        const double tangent_step = fabs(value_of_rec(delta.d_));
        const double tangent_scale = fabs(value_of_rec(f.d_ / f.val_));
        if (tangent_step > precision * std::max(1.0, tangent_scale)) {
          continue;
        }
      }
      break;
    }
  }
  return log_prefactor - log(f);
}

/**
 * Compute log(Q(a,z)) and d/da log(Q(a,z)) for z > a + 1 from the
 * continued fraction of log_q_gamma_cf, with the derivative with respect
 * to a carried through the recurrence. Both results stay in log space, so
 * they are finite where Q(a,z) underflows.
 *
 * @tparam T_a Type of shape parameter a (double or autodiff types)
 * @tparam T_z Type of value parameter z (double or autodiff types)
 * @param a Shape parameter
 * @param z Value at which to evaluate
 * @param precision Convergence threshold for the value and the derivative
 * @param max_steps Maximum number of continued fraction iterations
 * @return pair of log(Q(a,z)) and d/da log(Q(a,z))
 */
template <typename T_a, typename T_z>
inline std::pair<return_type_t<T_a, T_z>, return_type_t<T_a, T_z>>
log_q_gamma_cf_and_da(const T_a& a, const T_z& z,
                      double precision = LOG_Q_GAMMA_CF_PRECISION,
                      int max_steps = 250) {
  using T_return = return_type_t<T_a, T_z>;
  const T_return log_z = log(z);
  const T_return log_prefactor = a * log_z - z - lgamma(a);
  const T_return d_log_prefactor = log_z - digamma(a);

  // Each quantity x is carried with its derivative dx/da.
  const T_return b_init = z + 1.0 - a;
  T_return C = b_init;
  T_return dC = -1.0;
  if (fabs(value_of_rec(C)) < EPSILON) {
    C = EPSILON;
    dC = 0.0;
  }
  T_return D = 0.0;
  T_return dD = 0.0;
  T_return f = C;
  T_return df = dC;
  for (int i = 1; i <= max_steps; ++i) {
    const T_a an = -i * (i - a);
    const T_return b = b_init + 2.0 * i;
    dD = -1.0 + i * D + an * dD;
    D = b + an * D;
    if (fabs(value_of_rec(D)) < EPSILON) {
      D = EPSILON;
      dD = 0.0;
    }
    dC = -1.0 + (i * C - an * dC) / (C * C);
    C = b + an / C;
    if (fabs(value_of_rec(C)) < EPSILON) {
      C = EPSILON;
      dC = 0.0;
    }
    dD = -dD / (D * D);
    D = inv(D);
    const T_return delta = C * D;
    const T_return d_delta = dC * D + C * dD;
    df = df * delta + f * d_delta;
    f *= delta;
    const double tangent_scale = fabs(value_of_rec(df / f));
    if (fabs(value_of_rec(delta) - 1.0) < precision
        && fabs(value_of_rec(d_delta))
               <= precision * std::max(1.0, tangent_scale)) {
      break;
    }
  }
  return {log_prefactor - log(f), d_log_prefactor - df / f};
}

/**
 * Lower limit of log(Q(a,z)) for d/da log(Q) = (dQ/da) / Q. Below it Q or
 * dQ/da can be subnormal or zero, and log_q_gamma_cf_and_da must be used.
 */
inline double log_q_gamma_ratio_min() {
  static const double log_min
      = std::log(std::numeric_limits<double>::min()) - LOG_EPSILON;
  return log_min;
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
 * @param precision convergence threshold, default of sqrt(machine_epsilon)
 * @param max_steps maximum iterations for continued fraction
 * @return structure containing log(Q(a,z)) and d/da log(Q(a,z))
 */
template <typename T_a, typename T_z>
inline std::pair<return_type_t<T_a, T_z>, return_type_t<T_a, T_z>>
log_gamma_q_dgamma(const T_a& a, const T_z& z,
                   double precision = internal::LOG_Q_GAMMA_CF_PRECISION,
                   int max_steps = 250) {
  using T_return = return_type_t<T_a, T_z>;
  const double a_val = value_of(a);
  const double z_val = value_of(z);
  // For z > a + 1, use continued fraction for better numerical stability
  if (z_val > a_val + 1.0) {
    std::pair<T_return, T_return> result{
        internal::log_q_gamma_cf(a_val, z_val, precision, max_steps), 0.0};
    if (value_of_rec(result.first) < internal::log_q_gamma_ratio_min()) {
      const std::pair<double, double> log_q_and_da
          = internal::log_q_gamma_cf_and_da(a_val, z_val, precision, max_steps);
      result.first = log_q_and_da.first;
      result.second = log_q_and_da.second;
      return result;
    }
    // For gradient, use: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const T_return Q_val = exp(result.first);
    const double dQ_da
        = grad_reg_inc_gamma(a_val, z_val, tgamma(a_val), digamma(a_val));
    result.second = dQ_da / Q_val;
    return result;
  } else {
    // For z <= a + 1, use log1m(P(a,z)) for better numerical accuracy
    const double P_val = gamma_p(a_val, z_val);
    std::pair<T_return, T_return> result{log1m(P_val), 0.0};
    // Gradient: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const T_return Q_val = exp(result.first);
    if (Q_val > 0) {
      const double dQ_da
          = grad_reg_inc_gamma(a_val, z_val, tgamma(a_val), digamma(a_val));
      result.second = dQ_da / Q_val;
    } else {
      // Fallback if Q rounds to zero - use asymptotic approximation
      result.second = log(z_val) - digamma(a_val);
    }
    return result;
  }
}

}  // namespace math
}  // namespace stan

#endif
