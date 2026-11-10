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
#include <stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
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
 * Compute log(P(a,z)) and d/da log(P(a,z)) for z < a + 1 from the power
 * series P(a,z) = z^a exp(-z) / Gamma(a + 1) * sum_k t_k, with t_0 = 1 and
 * t_k = t_{k-1} z / (a + k), in log space. The derivative of log(t_k) with
 * respect to a is carried with every term. The series stops when a bound
 * on the remaining terms is below machine epsilon for the sum and for its
 * derivative.
 *
 * @tparam T_a Type of shape parameter a (double or autodiff types)
 * @tparam T_z Type of value parameter z (double or autodiff types)
 * @param a Shape parameter
 * @param z Value at which to evaluate, 0 < z < a + 1
 * @param max_steps Maximum number of terms
 * @return pair of log(P(a,z)) and d/da log(P(a,z))
 */
template <typename T_a, typename T_z>
inline std::pair<return_type_t<T_a, T_z>, return_type_t<T_a, T_z>>
log_p_gamma_series_and_da(const T_a& a, const T_z& z, int max_steps = 100000) {
  using T_return = return_type_t<T_a, T_z>;
  const double a_val = value_of_rec(a);
  const double z_val = value_of_rec(z);
  T_return t = 1.0;
  T_return d_log_t = 0.0;
  T_return sum = 1.0;
  T_return d_sum = 0.0;
  for (int k = 1; k <= max_steps; ++k) {
    t *= z / (a + k);
    d_log_t -= inv(a + k);
    sum += t;
    d_sum += t * d_log_t;
    // the remaining terms sum to at most t q / (1 - q)
    const double q = z_val / (a_val + k + 1);
    if (q < 1) {
      const double tail = value_of_rec(t) * q / (1 - q);
      if (tail <= EPSILON * value_of_rec(sum)
          && tail * (fabs(value_of_rec(d_log_t)) + 1)
                 <= EPSILON * fabs(value_of_rec(d_sum))) {
        break;
      }
    }
  }
  const T_return log_z = log(z);
  return {a * log_z - z - lgamma(a + 1.0) + log(sum),
          log_z - digamma(a + 1.0) + d_sum / sum};
}

/**
 * Result of log_gamma_tail: a log probability of the regularized incomplete
 * gamma function and its derivatives with respect to a and z.
 */
template <typename T>
struct log_gamma_tail_result {
  T log_prob;
  T d_a;
  T d_z;
};

/**
 * Compute log(P(a,z)) (lower tail) or log(Q(a,z)) (upper tail) with its
 * derivatives with respect to a and z. Every result stays in log space, so
 * it is finite where the probability underflows.
 *
 * For z > a + 1, log(Q) and its derivative come from the continued
 * fraction, and log(P) = log1m_exp(log(Q)). Otherwise P comes from gamma_p;
 * where P is below DBL_MIN / EPSILON, log(P) and its derivative come from
 * the power series in log space.
 *
 * @tparam Upper true for log(Q), false for log(P)
 * @tparam NeedDa compute the derivative with respect to a
 * @tparam NeedDz compute the derivative with respect to z
 * @tparam T_a Type of shape parameter a (double or autodiff types)
 * @tparam T_z Type of value parameter z (double or autodiff types)
 * @param a Shape parameter, a > 0
 * @param z Value at which to evaluate, 0 <= z < inf
 * @return log probability, d/da and d/dz of it (0 where not computed)
 */
template <bool Upper, bool NeedDa, bool NeedDz, typename T_a, typename T_z>
inline log_gamma_tail_result<return_type_t<T_a, T_z>> log_gamma_tail(
    const T_a& a, const T_z& z) {
  using T_return = return_type_t<T_a, T_z>;
  const double a_val = value_of_rec(a);
  const double z_val = value_of_rec(z);
  log_gamma_tail_result<T_return> out{0.0, 0.0, 0.0};
  if (z_val == 0) {
    if constexpr (Upper) {
      // Q(a, 0) = 1; d/dz log(Q) at 0 is minus the density at 0
      if constexpr (NeedDz) {
        out.d_z = a_val > 1 ? 0.0 : (a_val == 1 ? -1.0 : NEGATIVE_INFTY);
      }
    } else {
      out.log_prob = NEGATIVE_INFTY;
    }
    return out;
  }
  if (z_val == INFTY) {
    // P(a, inf) = 1, Q(a, inf) = 0; the derivatives are treated as zero
    if constexpr (Upper) {
      out.log_prob = NEGATIVE_INFTY;
    }
    return out;
  }
  if (z_val > a_val + 1.0) {
    T_return log_q;
    T_return d_a_log_q = 0.0;
    if constexpr (NeedDa) {
      std::pair<T_return, T_return> cf = log_q_gamma_cf_and_da(a, z);
      log_q = cf.first;
      d_a_log_q = cf.second;
    } else {
      log_q = log_q_gamma_cf(a, z);
    }
    if constexpr (Upper) {
      out.log_prob = log_q;
      out.d_a = d_a_log_q;
    } else {
      out.log_prob = log1m_exp(log_q);
      if constexpr (NeedDa) {
        out.d_a = -exp(log_q - out.log_prob) * d_a_log_q;
      }
    }
  } else {
    const T_return p = gamma_p(a, z);
    if constexpr (Upper) {
      out.log_prob = log1m(p);
      if constexpr (NeedDa) {
        out.d_a = -grad_reg_lower_inc_gamma(a, z) / (1.0 - p);
      }
    } else if (value_of_rec(p)
               >= std::numeric_limits<double>::min() / EPSILON) {
      out.log_prob = log(p);
      if constexpr (NeedDa) {
        out.d_a = grad_reg_lower_inc_gamma(a, z) / p;
      }
    } else {
      std::pair<T_return, T_return> series = log_p_gamma_series_and_da(a, z);
      out.log_prob = series.first;
      out.d_a = series.second;
    }
  }
  if constexpr (NeedDz) {
    // the density of the gamma distribution at z, divided by P or Q
    const T_return log_density = (a - 1.0) * log(z) - z - lgamma(a);
    const T_return ratio = exp(log_density - out.log_prob);
    if constexpr (Upper) {
      out.d_z = -ratio;
    } else {
      out.d_z = ratio;
    }
  }
  return out;
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
  // For z > a + 1, use the continued fraction, with the derivative carried
  // through it in log space; it stays finite where Q underflows
  if (z_val > a_val + 1.0) {
    const std::pair<double, double> log_q_and_da
        = internal::log_q_gamma_cf_and_da(a_val, z_val, precision, max_steps);
    return {log_q_and_da.first, log_q_and_da.second};
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
