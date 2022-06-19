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
#include <cmath>

namespace stan {
namespace math {
namespace internal {
/**
 * Returns the intermediate calculations for the gradient of the
 * hypergeometric function (2F1).
 *
 * @tparam T scalar type or `var_value`
 * @tparam T_StdVec std::array of doubles or `var_values`
 * @tparam T_StdVecInt array of int
 * @tparam T_int integer
 * @param[out] log_g_old log of previous gradient values
 * @param[out] log_g_old_sign sign of previous gradient values
 * @param[out] log_t_new log of power series calculation
 * @param[out] log_t_new_sign sign of power series calcuation
 * @param[in] p intermediate power series value
 * @param[in] log_z see hypergeometric function 2F1 definition
 * @param[in] log_t_old log of previous power series calculation
 * @param[in] log_t_old_sign log of previous power series calculation
 * @param[in] k loop iteration
 * @param[in] a1 see hypergeometric function 2F1 definition
 * @param[in] a2 see hypergeometric function 2F1 definition
 * @param[in] b1 see hypergeometric function 2F1 definition
 */
template <typename T, typename T_StdVec, typename T_StdVecInt, typename T_int>
void calc_lambda(T_StdVec& log_g_old, T_StdVecInt& log_g_old_sign, T& log_t_new,
                 T_int& log_t_new_sign, const T& p, const T& log_z,
                 const T& log_t_old, const T_int& log_t_old_sign, const int k,
                 const T& a1, const T& a2, const T& b1) {
  using ret_t = return_type_t<T, T_int, T_StdVec, T_StdVecInt>;

  log_t_new += log(fabs(p)) + log_z;
  log_t_new_sign = sign(value_of_rec(p)) * log_t_new_sign;

  Eigen::Array<ret_t, Eigen::Dynamic, 1> hyper_args(3);
  hyper_args << a1 + k, a2 + k, -(b1 + k);

  Eigen::Array<ret_t, Eigen::Dynamic, 1> term
      = log_g_old_sign * log_t_old_sign * exp(log_g_old - log_t_old)
        + inv(hyper_args);

  log_g_old = log_t_new + log(abs(term));
  log_g_old_sign = sign(value_of_rec(term)) * log_t_new_sign;
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
 * @tparam T1 scalar or `var_value`
 * @tparam T2 scalar or `var_value`
 * @tparam T3 scalar or `var_value`
 * @tparam T_z scalar type
 * @param[out] g_a1 reference to gradient of 2F1 w.r.t. a1, result
 * @param[out] g_a2 reference to gradient of 2F1 w.r.t. a2, result
 * @param[out] g_b1 reference to gradient of 2F1 w.r.t. b1, result
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 */
template <typename T1, typename T2, typename T3, typename T_z>
void grad_2F1(T1& g_a1, T2& g_a2, T3& g_b1, const T1& a1, const T2& a2,
              const T3& b1, const T_z& z, double precision = 1e-14,
              int max_steps = 1e6) {
  check_2F1_converges("grad_2F1", a1, a2, b1, z);

  using stan::math::value_of_rec;
  using std::max;
  using ret_t = return_type_t<T1, T2, T3, T_z>;

  if (z == 0) {
    return;
  }

  Eigen::Array<ret_t, Eigen::Dynamic, 1> g(3);
  g << 0.0, 0.0, 0.0;

  Eigen::Array<ret_t, Eigen::Dynamic, 1> log_g_old(3);
  log_g_old << NEGATIVE_INFTY, NEGATIVE_INFTY, NEGATIVE_INFTY;

  ret_t log_t_old = 0.0;
  ret_t log_t_new = 0.0;
  int sign_z = sign(z);
  ret_t log_z = log(abs(z));

  ret_t log_precision = log(precision);
  ret_t log_t_new_sign = 1.0;
  ret_t log_t_old_sign = 1.0;

  Eigen::Array<ret_t, Eigen::Dynamic, 1> log_g_old_sign(3);
  log_g_old_sign << 1., 1., 1.;

  int sign_zk = sign_z;

  for (int k = 0; k <= max_steps; ++k) {
    ret_t p = ((a1 + k) * (a2 + k) / ((b1 + k) * (1 + k)));
    if (p == 0) {
      return;
    }
    internal::calc_lambda(log_g_old, log_g_old_sign, log_t_new, log_t_new_sign,
                          p, log_z, log_t_old, log_t_old_sign, k, a1, a2, b1);

    g += log_g_old_sign * exp(log_g_old) * sign_zk;
    g_a1 = g(0);
    g_a2 = g(1);
    g_b1 = g(2);

    if (log_g_old(0)
            <= max(log(abs(value_of_rec(g(0)))) + log_precision, log_precision)
        && log_g_old(1) <= max(
               log(std::abs(value_of_rec(g(1)))) + log_precision, log_precision)
        && log_g_old(2) <= max(log(abs(value_of_rec(g(2)))) + log_precision,
                               log_precision)) {
      return;
    }
    log_t_old = log_t_new;
    log_t_old_sign = log_t_new_sign;
    sign_zk *= sign_z;
  }
  throw_domain_error("grad_2F1", "k (internal counter)", max_steps, "exceeded ",
                     " iterations, hypergeometric function gradient "
                     "did not converge.");
  return;
}

}  // namespace math
}  // namespace stan
#endif
