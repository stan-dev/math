#ifndef STAN_MATH_PRIM_FUN_GRAD_2F1_HPP
#define STAN_MATH_PRIM_FUN_GRAD_2F1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {
template <typename T_p, typename T_log_z, typename T_input, typename T_log,
          typename T_array, typename T_array_int, typename T_sign>
void calc_lambda(const T_p& p, const T_log_z& log_z, T_array& log_g_old,
                 T_array_int& log_g_old_sign, T_log& log_t_old,
                 T_sign& log_t_old_sign, T_log& log_t_new,
                 T_sign& log_t_new_sign, int k, T_input& a1, T_input& a2,
                 T_input& b1) {
  using T_plain = return_type_t<T_log, T_sign, T_array, T_array_int>;

  log_t_new += log(fabs(p)) + log_z;
  log_t_new_sign = p >= 0.0 ? log_t_new_sign : -log_t_new_sign;

  T_plain term
      = log_g_old_sign[0] * log_t_old_sign * exp(log_g_old[0] - log_t_old)
        + inv(a1 + k);
  log_g_old[0] = log_t_new + log(fabs(term));
  log_g_old_sign[0] = term >= 0.0 ? log_t_new_sign : -log_t_new_sign;

  term = log_g_old_sign[1] * log_t_old_sign * exp(log_g_old[1] - log_t_old)
         + inv(a2 + k);
  log_g_old[1] = log_t_new + log(fabs(term));
  log_g_old_sign[1] = term >= 0.0 ? log_t_new_sign : -log_t_new_sign;

  term = log_g_old_sign[2] * log_t_old_sign * exp(log_g_old[2] - log_t_old)
         - inv(b1 + k);
  log_g_old[2] = log_t_new + log(fabs(term));
  log_g_old_sign[2] = term >= 0.0 ? log_t_new_sign : -log_t_new_sign;
}
}  // namespace internal

/**
 * Gradients of the hypergeometric function, 2F1.
 *
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * This power-series representation converges for all gradients
 * under the same conditions as the 2F1 function itself.
 *
 * @tparam T type of arguments and result
 * @param[out] g_a1 g_a1 reference to gradient of 2F1 w.r.t. a1, result.
 * @param[out] g_a2 g_a2 reference to gradient of 2F1 w.r.t. a2, result.
 * @param[out] g_b1 g_b1 reference to gradient of 2F1 w.r.t. b1, result.
 * @param[in] a1 a1 see generalized hypergeometric function definition.
 * @param[in] a2 a2 see generalized hypergeometric function definition.
 * @param[in] b1 b1 see generalized hypergeometric function definition.
 * @param[in] z z see generalized hypergeometric function definition.
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum at.
 * @param[in] max_steps number of steps to take.
 */
template <typename T1, typename T2, typename T3, typename T_z>
void grad_2F1(T1& g_a1, T2& g_a2, T3& g_b1, const T1& a1, const T2& a2,
              const T3& b1, const T_z& z, double precision = 1e-14,
              int max_steps = 1e6) {
  check_2F1_converges("grad_2F1", a1, a2, b1, z);

  using stan::math::value_of_rec;
  using std::exp;
  using std::fabs;
  using std::log;
  using std::max;
  using TP = return_type_t<T1, T2, T3, T_z>;

  g_a1 = 0.0;
  g_a2 = 0.0;
  g_b1 = 0.0;

  TP log_g_old[3];
  for (auto& i : log_g_old) {
    i = NEGATIVE_INFTY;
  }

  TP log_t_old = 0.0;
  TP log_t_new = 0.0;
  int sign_z = sign(z);
  TP log_z = log(fabs(z));

  TP log_precision = log(precision);
  TP log_t_new_sign = 1.0;
  TP log_t_old_sign = 1.0;
  TP log_g_old_sign[3];
  for (TP& x : log_g_old_sign) {
    x = 1.0;
  }

  int sign_zk = 1;

  if (sign_z > 0) {
    for (int k = 0; k <= max_steps; ++k) {
      TP p = ((a1 + k) * (a2 + k) / ((b1 + k) * (1 + k)));
      if (p == 0) {
        return;
      }

      internal::calc_lambda(p, log_z, log_g_old, log_g_old_sign, log_t_old,
                            log_t_old_sign, log_t_new, log_t_new_sign, k, a1,
                            a2, b1);

      g_a1 += log_g_old_sign[0] > 0 ? exp(log_g_old[0]) : -exp(log_g_old[0]);
      g_a2 += log_g_old_sign[1] > 0 ? exp(log_g_old[1]) : -exp(log_g_old[1]);
      g_b1 += log_g_old_sign[2] > 0 ? exp(log_g_old[2]) : -exp(log_g_old[2]);

      if (log_g_old[0] <= std::max(
              std::log(std::abs(value_of_rec(g_a1))) + log_precision,
              log_precision)
          && log_g_old[1] <= std::max(
                 std::log(std::abs(value_of_rec(g_a2))) + log_precision,
                 log_precision)
          && log_g_old[2] <= std::max(
                 std::log(std::abs(value_of_rec(g_b1))) + log_precision,
                 log_precision)) {
        return;
      }

      log_t_old = log_t_new;
      log_t_old_sign = log_t_new_sign;
    }
  } else {
    for (int k = 0; k <= max_steps; ++k) {
      TP p = ((a1 + k) * (a2 + k) / ((b1 + k) * (1 + k)));
      if (p == 0) {
        return;
      }
      internal::calc_lambda(p, log_z, log_g_old, log_g_old_sign, log_t_old,
                            log_t_old_sign, log_t_new, log_t_new_sign, k, a1,
                            a2, b1);

      if (sign_zk > 0) {
        log_t_new += log(fabs(p)) + log_z;
        log_t_new_sign = p >= 0.0 ? log_t_new_sign : -log_t_new_sign;

        g_a1 += log_g_old_sign[0] > 0 ? -exp(log_g_old[0]) : exp(log_g_old[0]);
        g_a2 += log_g_old_sign[1] > 0 ? -exp(log_g_old[1]) : exp(log_g_old[1]);
        g_b1 += log_g_old_sign[2] > 0 ? -exp(log_g_old[2]) : exp(log_g_old[2]);
      } else {
        g_a1 += log_g_old_sign[0] > 0 ? exp(log_g_old[0]) : -exp(log_g_old[0]);
        g_a2 += log_g_old_sign[1] > 0 ? exp(log_g_old[1]) : -exp(log_g_old[1]);
        g_b1 += log_g_old_sign[2] > 0 ? exp(log_g_old[2]) : -exp(log_g_old[2]);
      }
      if (log_g_old[0] <= std::max(
              std::log(std::abs(value_of_rec(g_a1))) + log_precision,
              log_precision)
          && log_g_old[1] <= std::max(
                 std::log(std::abs(value_of_rec(g_a2))) + log_precision,
                 log_precision)
          && log_g_old[2] <= std::max(
                 std::log(std::abs(value_of_rec(g_b1))) + log_precision,
                 log_precision)) {
        return;
      }
      log_t_old = log_t_new;
      log_t_old_sign = log_t_new_sign;
      sign_zk *= sign_z;
    }
  }
  throw_domain_error("grad_2F1", "k (internal counter)", max_steps, "exceeded ",
                     " iterations, hypergeometric function gradient "
                     "did not converge.");
  return;
}

}  // namespace math
}  // namespace stan
#endif
