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
#include <cmath>

namespace stan {
namespace math {

namespace internal {

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
                                              double precision = 1.49012e-08,
                                              int max_steps = 250) {
  using T_return = return_type_t<T_a, T_z>;
  const auto log_prefactor = a * log(z) - z - lgamma(a);

  auto b_init = z + 1.0 - a;
  auto C = (fabs(value_of(b_init)) >= EPSILON) ? b_init : std::decay_t<decltype(b_init)>(EPSILON);
  auto D = 0.0;
  auto f = C;
  for (int i = 1; i <= max_steps; ++i) {
    T_a an = -i * (i - a);
    const auto b = b_init + 2.0 * i;
    D = b + an * D;
    D = (fabs(value_of(D)) >= EPSILON) ? D : std::decay_t<decltype(D)>(EPSILON);
    C = b + an / C;
    C = (fabs(value_of(C)) >= EPSILON) ? C : std::decay_t<decltype(C)>(EPSILON);
    D = inv(D);
    auto delta = C * D;
    f *= delta;
    const auto delta_m1 = fabs(value_of(delta) - 1.0);
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
 * @param precision convergence threshold, default of sqrt(machine_epsilon)
 * @param max_steps maximum iterations for continued fraction
 * @return structure containing log(Q(a,z)) and d/da log(Q(a,z))
 */
template <typename T_a, typename T_z>
inline auto log_gamma_q_dgamma(
    const T_a& a, const T_z& z, double precision = 1.49012e-08, int max_steps = 250) {
  using T_return = return_type_t<T_a, T_z>;
  const auto a_val = value_of(a);
  const auto z_val = value_of(z);
  // For z > a + 1, use continued fraction for better numerical stability
  if (z_val > a_val + 1.0) {
    std::pair<T_return, T_return> result{internal::log_q_gamma_cf(a_val, z_val, precision, max_steps), 0.0};
    // For gradient, use: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const auto Q_val = exp(result.first);
    const auto dQ_da
        = grad_reg_inc_gamma(a_val, z_val, tgamma(a_val), digamma(a_val));
    result.second = dQ_da / Q_val;
    return result;
  } else {
    // For z <= a + 1, use log1m(P(a,z)) for better numerical accuracy
    const auto P_val = gamma_p(a_val, z_val);
    std::pair<T_return, T_return> result{log1m(P_val), 0.0};
    // Gradient: d/da log(Q) = (1/Q) * dQ/da
    // grad_reg_inc_gamma computes dQ/da
    const auto Q_val = exp(result.first);
    if (Q_val > 0) {
      const auto dQ_da
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
