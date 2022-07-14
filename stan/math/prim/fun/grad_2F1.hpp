#ifndef STAN_MATH_PRIM_FUN_GRAD_2F1_HPP
#define STAN_MATH_PRIM_FUN_GRAD_2F1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/hypergeometric_2F1.hpp>
#include <cmath>

namespace stan {
namespace math {

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
 * @tparam T1 scalar or `var_value`
 * @tparam T2 scalar or `var_value`
 * @tparam T3 scalar or `var_value`
 * @tparam T_z scalar type
 * @param[out] g_a1 reference to gradient of 2F1 w.r.t. a1, result
 * @param[out] g_a2 reference to gradient of 2F1 w.r.t. a2, result
 * @param[out] g_b1 reference to gradient of 2F1 w.r.t. b1, result
 * @param[out] g_z reference to gradient of 2F1 w.r.t. z, result
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 */
template <bool calc_a1, bool calc_a2, bool calc_b1, bool calc_z,
          typename T1, typename T2, typename T3, typename T_z>
void grad_2F1(T1& g_a1, T2& g_a2, T3& g_b1, T_z& g_z, const T1& a1,
              const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  check_2F1_converges("grad_2F1", a1, a2, b1, z);

  using stan::math::value_of_rec;
  using std::max;
  using ret_t = return_type_t<T1, T2, T3, T_z>;

  if (calc_z) {
    auto hyper_2f1_dz = hypergeometric_2F1(a1 + 1.0, a2 + 1.0, b1 + 1.0, z);
    g_z = (a1 * a2 * hyper_2f1_dz) / b1;
  }

  if (z == 0) {
    return;
  }

  if (!(calc_a1 || calc_a2 || calc_b1)) {
    return;
  }

  Eigen::Array<ret_t, Eigen::Dynamic, 1> g(3);
  g << 0.0, 0.0, 0.0;


  Eigen::Array<ret_t, Eigen::Dynamic, 1> log_g_old(3);
  log_g_old << NEGATIVE_INFTY, NEGATIVE_INFTY, NEGATIVE_INFTY;

  ret_t log_t_old = 0.0;
  ret_t log_t_new = 0.0;
  int sign_z = sign(z);
  auto log_z = log(abs(z));

  double log_precision = log(precision);
  int log_t_new_sign = 1.0;
  int log_t_old_sign = 1.0;

  Eigen::Array<ret_t, Eigen::Dynamic, 1> log_g_old_sign(3);
  log_g_old_sign << 1., 1., 1.;

  int sign_zk = sign_z;
  int k = 0;
  ret_t inner_diff = 1;
  Eigen::Matrix<ret_t, -1, 1> g_current(3);
  g_current << 0, 0, 0;

  while (inner_diff > precision && k < max_steps) {
    ret_t p = ((a1 + k) * (a2 + k) / ((b1 + k) * (1 + k)));
    if (p == 0) {
      return;
    }
    log_t_new += log(fabs(p)) + log_z;
    log_t_new_sign = sign(value_of_rec(p)) * log_t_new_sign;

    if (calc_a1) {
      ret_t term_a1 = log_g_old_sign(0)
                      * log_t_old_sign
                      * exp(log_g_old(0) - log_t_old)
                      + inv(a1 + k);
      log_g_old(0) = log_t_new + log(abs(term_a1));
      log_g_old_sign(0) = sign(value_of_rec(term_a1)) * log_t_new_sign;
      g_current(0) = log_g_old_sign(0) * exp(log_g_old(0)) * sign_zk;
      g(0) += g_current(0);
    }

    if (calc_a2) {
      ret_t term_a2 = log_g_old_sign(1)
                      * log_t_old_sign
                      * exp(log_g_old(1) - log_t_old)
                      + inv(a2 + k);
      log_g_old(1) = log_t_new + log(abs(term_a2));
      log_g_old_sign(1) = sign(value_of_rec(term_a2)) * log_t_new_sign;
      g_current(1) = log_g_old_sign(1) * exp(log_g_old(1)) * sign_zk;
      g(1) += g_current(1);
    }

    if (calc_b1) {
      ret_t term_b1 = log_g_old_sign(2)
                      * log_t_old_sign
                      * exp(log_g_old(2) - log_t_old)
                      + inv(-(b1 + k));
      log_g_old(2) = log_t_new + log(abs(term_b1));
      log_g_old_sign(2) = sign(value_of_rec(term_b1)) * log_t_new_sign;
      g_current(2) = log_g_old_sign(2) * exp(log_g_old(2)) * sign_zk;
      g(2) += g_current(2);
    }

    inner_diff = g_current.maxCoeff();

    log_t_old = log_t_new;
    log_t_old_sign = log_t_new_sign;
    sign_zk *= sign_z;
    ++k;
  }

  g_a1 = g(0);
  g_a2 = g(1);
  g_b1 = g(2);

  if (k > max_steps) {
    throw_domain_error("grad_2F1", "k (internal counter)", max_steps, "exceeded ",
                      " iterations, hypergeometric function gradient "
                      "did not converge.");
  }
  return;
}

template <typename T1, typename T2, typename T3, typename T_z>
void grad_2F1(T1& g_a1, T2& g_a2, T3& g_b1, T_z& g_z, const T1& a1,
              const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  grad_2F1<true, true, true, true>(g_a1, g_a2, g_b1, g_z, a1, a2, b1, z, precision, max_steps);
}

}  // namespace math
}  // namespace stan
#endif
