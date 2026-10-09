#ifndef STAN_MATH_PRIM_FUN_GRAD_F32_HPP
#define STAN_MATH_PRIM_FUN_GRAD_F32_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/hypergeometric_3F2_tail_bound.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_nonpositive_integer.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <array>
#include <cmath>

namespace stan {
namespace math {

/**
 * Gradients of the hypergeometric function, 3F2.
 *
 * Calculate the gradients of the hypergeometric function (3F2)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * A terminating series (a1, a2 or a3 a non-positive integer) is a
 * polynomial. Its sum stops at the first zero term, or earlier when bounds
 * of the sums of all remaining terms of the series and of every requested
 * gradient are below 1e-17 times the absolute values of their partial sums;
 * <code>precision</code> does not apply to it.
 *
 * This power-series representation converges for all gradients
 * under the same conditions as the 3F2 function itself.
 *
 * @tparam grad_a1 boolean indicating if gradient with respect to a1 is required
 * @tparam grad_a2 boolean indicating if gradient with respect to a2 is required
 * @tparam grad_a3 boolean indicating if gradient with respect to a3 is required
 * @tparam grad_b1 boolean indicating if gradient with respect to b1 is required
 * @tparam grad_b2 boolean indicating if gradient with respect to b2 is required
 * @tparam grad_z boolean indicating if gradient with respect to z is required
 * @tparam T1 a scalar type
 * @tparam T2 a scalar type
 * @tparam T3 a scalar type
 * @tparam T4 a scalar type
 * @tparam T5 a scalar type
 * @tparam T6 a scalar type
 * @tparam T7 a scalar type
 * @tparam T8 a scalar type
 * @param[out] g g pointer to array of six values of type T, result.
 * @param[in] a1 a1 see generalized hypergeometric function definition.
 * @param[in] a2 a2 see generalized hypergeometric function definition.
 * @param[in] a3 a3 see generalized hypergeometric function definition.
 * @param[in] b1 b1 see generalized hypergeometric function definition.
 * @param[in] b2 b2 see generalized hypergeometric function definition.
 * @param[in] z z see generalized hypergeometric function definition.
 * @param[in] precision precision of the infinite sum: the sum stops after a
 *   term below <code>precision</code> in absolute value. Not used for a
 *   terminating series.
 * @param[in] max_steps number of steps to take
 */
template <bool grad_a1 = true, bool grad_a2 = true, bool grad_a3 = true,
          bool grad_b1 = true, bool grad_b2 = true, bool grad_z = true,
          typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8 = double>
inline void grad_F32(T1* g, const T2& a1, const T3& a2, const T4& a3,
                     const T5& b1, const T6& b2, const T7& z,
                     const T8& precision = 1e-6, int max_steps = 1e5) {
  check_3F2_converges("grad_F32", a1, a2, a3, b1, b2, z);

  for (int i = 0; i < 6; ++i) {
    g[i] = 0.0;
  }

  T1 log_g_old[6];
  for (auto& x : log_g_old) {
    x = NEGATIVE_INFTY;
  }

  T1 log_t_old = 0.0;
  T1 log_t_new = 0.0;

  T7 log_z = log(z);

  T1 log_t_new_sign = 1.0;
  T1 log_t_old_sign = 1.0;
  T1 log_g_old_sign[6];
  for (int i = 0; i < 6; ++i) {
    log_g_old_sign[i] = 1.0;
  }
  // A terminating series stops when bounds of the sums of all remaining
  // terms of the series and of every requested gradient are negligible
  // against the partial sums.
  const bool is_polynomial = is_nonpositive_integer(a1)
                             || is_nonpositive_integer(a2)
                             || is_nonpositive_integer(a3);
  constexpr double polynomial_precision = 1e-17;
  const std::array<double, 3> a_val{value_of_rec(a1), value_of_rec(a2),
                                    value_of_rec(a3)};
  const std::array<double, 2> b_val{value_of_rec(b1), value_of_rec(b2)};
  const double abs_z = std::fabs(value_of_rec(z));
  // index of the last term that can be nonzero
  double last = max_steps + 1.0;
  for (double a_i : a_val) {
    if (a_i <= 0.0 && a_i == std::floor(a_i)) {
      last = std::fmin(last, -a_i);
    }
  }
  auto is_negligible = [](double increment, const auto& sum) {
    return increment <= polynomial_precision * std::fabs(value_of_rec(sum));
  };
  double f_sum = 1.0;
  std::array<T1, 6> term{0};
  for (int k = 0; k <= max_steps; ++k) {
    // A numerator parameter that has reached zero ends the series (a
    // polynomial). Stop before the ratio, which is 0 / 0 when a denominator
    // parameter reaches zero at the same k.
    if (a1 + k == 0 || a2 + k == 0 || a3 + k == 0) {
      return;
    }
    T1 p = (a1 + k) * (a2 + k) * (a3 + k) / ((b1 + k) * (b2 + k) * (1 + k));
    if (p == 0) {
      return;
    }

    log_t_new += log(fabs(p)) + log_z;
    log_t_new_sign = p >= 0.0 ? log_t_new_sign : -log_t_new_sign;
    bool negligible = false;
    double t_new = 0.0;
    if (is_polynomial) {
      t_new = std::exp(value_of_rec(log_t_new));
      f_sum += value_of_rec(log_t_new_sign) * t_new;
      negligible = is_negligible(t_new, f_sum);
    }
    if constexpr (grad_a1) {
      term[0]
          = log_g_old_sign[0] * log_t_old_sign * exp(log_g_old[0] - log_t_old)
            + inv(a1 + k);
      log_g_old[0] = log_t_new + log(fabs(term[0]));
      log_g_old_sign[0] = term[0] >= 0.0 ? log_t_new_sign : -log_t_new_sign;
      const T1 increment = exp(log_g_old[0]);
      g[0] += log_g_old_sign[0] * increment;
      negligible = negligible && is_negligible(value_of_rec(increment), g[0]);
    }

    if constexpr (grad_a2) {
      term[1]
          = log_g_old_sign[1] * log_t_old_sign * exp(log_g_old[1] - log_t_old)
            + inv(a2 + k);
      log_g_old[1] = log_t_new + log(fabs(term[1]));
      log_g_old_sign[1] = term[1] >= 0.0 ? log_t_new_sign : -log_t_new_sign;
      const T1 increment = exp(log_g_old[1]);
      g[1] += log_g_old_sign[1] * increment;
      negligible = negligible && is_negligible(value_of_rec(increment), g[1]);
    }

    if constexpr (grad_a3) {
      term[2]
          = log_g_old_sign[2] * log_t_old_sign * exp(log_g_old[2] - log_t_old)
            + inv(a3 + k);
      log_g_old[2] = log_t_new + log(fabs(term[2]));
      log_g_old_sign[2] = term[2] >= 0.0 ? log_t_new_sign : -log_t_new_sign;
      const T1 increment = exp(log_g_old[2]);
      g[2] += log_g_old_sign[2] * increment;
      negligible = negligible && is_negligible(value_of_rec(increment), g[2]);
    }

    if constexpr (grad_b1) {
      term[3]
          = log_g_old_sign[3] * log_t_old_sign * exp(log_g_old[3] - log_t_old)
            - inv(b1 + k);
      log_g_old[3] = log_t_new + log(fabs(term[3]));
      log_g_old_sign[3] = term[3] >= 0.0 ? log_t_new_sign : -log_t_new_sign;
      const T1 increment = exp(log_g_old[3]);
      g[3] += log_g_old_sign[3] * increment;
      negligible = negligible && is_negligible(value_of_rec(increment), g[3]);
    }

    if constexpr (grad_b2) {
      term[4]
          = log_g_old_sign[4] * log_t_old_sign * exp(log_g_old[4] - log_t_old)
            - inv(b2 + k);
      log_g_old[4] = log_t_new + log(fabs(term[4]));
      log_g_old_sign[4] = term[4] >= 0.0 ? log_t_new_sign : -log_t_new_sign;
      const T1 increment = exp(log_g_old[4]);
      g[4] += log_g_old_sign[4] * increment;
      negligible = negligible && is_negligible(value_of_rec(increment), g[4]);
    }

    if constexpr (grad_z) {
      term[5]
          = log_g_old_sign[5] * log_t_old_sign * exp(log_g_old[5] - log_t_old)
            + inv(z);
      log_g_old[5] = log_t_new + log(fabs(term[5]));
      log_g_old_sign[5] = term[5] >= 0.0 ? log_t_new_sign : -log_t_new_sign;
      const T1 increment = exp(log_g_old[5]);
      g[5] += log_g_old_sign[5] * increment;
      negligible = negligible && is_negligible(value_of_rec(increment), g[5]);
    }

    if (is_polynomial) {
      if (negligible) {
        // The remaining terms t_{k + 2}, ..., t_last have the ratios r_j for
        // j = k + 1, ..., last - 1. A gradient term is t_j s_j, where s_j
        // (term[i] for j = k + 1) grows by at most c per step.
        const double n_left = last - k - 1;
        if (n_left <= 0) {
          return;
        }
        const double lo = k + 1;
        const double hi = last - 1;
        const auto sums = internal::hypergeometric_3F2_tail_sums(
            internal::hypergeometric_3F2_ratio_bound(a_val, b_val, abs_z, lo,
                                                     hi),
            n_left);
        auto tail_is_negligible = [&](const auto& s, double c,
                                      const auto& sum) {
          return is_negligible(
              t_new
                  * (std::fabs(value_of_rec(s)) * sums.first + c * sums.second),
              sum);
        };
        auto growth = [&](double p) {
          return 1.0 / internal::hypergeometric_3F2_min_abs(p, lo, hi);
        };
        bool done = is_negligible(t_new * sums.first, f_sum);
        if constexpr (grad_a1) {
          done = done && tail_is_negligible(term[0], growth(a_val[0]), g[0]);
        }
        if constexpr (grad_a2) {
          done = done && tail_is_negligible(term[1], growth(a_val[1]), g[1]);
        }
        if constexpr (grad_a3) {
          done = done && tail_is_negligible(term[2], growth(a_val[2]), g[2]);
        }
        if constexpr (grad_b1) {
          done = done && tail_is_negligible(term[3], growth(b_val[0]), g[3]);
        }
        if constexpr (grad_b2) {
          done = done && tail_is_negligible(term[4], growth(b_val[1]), g[4]);
        }
        if constexpr (grad_z) {
          done = done && tail_is_negligible(term[5], 1.0 / abs_z, g[5]);
        }
        if (done) {
          return;
        }
      }
    } else if (log_t_new <= log(precision)) {
      return;  // implicit abs
    }

    log_t_old = log_t_new;
    log_t_old_sign = log_t_new_sign;
  }
  throw_domain_error("grad_F32", "k (internal counter)", max_steps, "exceeded ",
                     " iterations, hypergeometric function gradient "
                     "did not converge.");
  return;
}

}  // namespace math
}  // namespace stan
#endif
