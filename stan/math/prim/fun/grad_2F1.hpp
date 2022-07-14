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
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/hypergeometric_2F1.hpp>
#include <cmath>
#include <boost/optional.hpp>

namespace stan {
namespace math {

namespace internal {
/**
 * Implementation function to calculate the gradients of the hypergeometric
 * function, 2F1.
 *
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * This power-series representation converges for all gradients
 * under the same conditions as the 2F1 function itself. As with the
 * hypergeometric_2F1 function, if the parameters do not meet convergence
 * criteria then the gradients are calculated using Euler's transformation.
 *
 * @tparam calc_a1 boolean for whether to calculate gradients w.r.t a1
 * @tparam calc_a2 boolean for whether to calculate gradients w.r.t a2
 * @tparam calc_b1 boolean for whether to calculate gradients w.r.t b1
 * @tparam calc_z boolean for whether to calculate gradients w.r.t z
 * @tparam Tg1 type of outcome variable to store a1 gradient
 * @tparam Tg2 type of outcome variable to store a2 gradient
 * @tparam Tg3 type of outcome variable to store b1 gradient
 * @tparam T_gz type of outcome variable to store z gradient
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
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
template <bool calc_a1, bool calc_a2, bool calc_b1, bool calc_z, typename Tg1,
          typename Tg2, typename Tg3, typename T_gz, typename T1, typename T2,
          typename T3, typename T_z,
          typename ScalarT = return_type_t<T1, T2, T3, T_z>,
          typename TupleT = std::tuple<T1, T2, T3, T_z>,
          typename OptT = boost::optional<TupleT>>
void grad_2F1_impl(Tg1& g_a1, Tg2& g_a2, Tg3& g_b1, T_gz& g_z, const T1& a1,
                   const T2& a2, const T3& b1, const T_z& z,
                   double precision = 1e-14, int max_steps = 1e6) {
  bool euler_transform = false;
  try {
    check_2F1_converges("hypergeometric_2F1", a1, a2, b1, z);
  } catch (...) {
    // Apply Euler's hypergeometric transformation if function
    // will not converge with current arguments
    check_2F1_converges("hypergeometric_2F1", b1 - a1, a2, b1, z / (z - 1));
    euler_transform = true;
  }

  ScalarT a1_t = euler_transform ? a2 : a1;
  ScalarT a2_t = euler_transform ? b1 - a1 : a2;
  ScalarT z_t = euler_transform ? z / (z - 1) : z;

  using stan::math::value_of_rec;
  using std::max;
  using ret_t = return_type_t<T1, T2, T3, T_z>;

  if (calc_z) {
    if (euler_transform) {
      auto hyper1 = hypergeometric_2F1(a1_t, a2_t, b1, z_t);
      auto hyper2 = hypergeometric_2F1(1 + a2, 1 - a1 + b1, 1 + b1, z_t);
      auto pre_mult = a2 * pow(1 - z, -1 - a2);
      g_z = a2 * pow(1 - z, -1 - a2) * hyper1
            + (a2 * (b1 - a1) * pow(1 - z, -a2)
               * (inv(z - 1) - z / square(z - 1)) * hyper2)
                  / b1;
    } else {
      auto hyper_2f1_dz = hypergeometric_2F1(a1 + 1.0, a2 + 1.0, b1 + 1.0, z);
      g_z = (a1 * a2 * hyper_2f1_dz) / b1;
    }
  }

  if (z == 0) {
    return;
  }

  if (!(calc_a1 || calc_a2 || calc_b1)) {
    return;
  }

  bool calc_a1_euler = calc_a1;
  bool calc_a2_euler = calc_a2;

  if (euler_transform) {
    // 'a' gradients under Euler transform are constructed using the gradients
    // of both elements, so need to compute both if any are required
    if (calc_a1 || calc_a2) {
      calc_a1_euler = true;
      calc_a2_euler = true;
    }
    // 'b' gradients under Euler transform require gradients from 'a2'
    if (calc_b1) {
      calc_a2_euler = true;
    }
  }

  Eigen::Array<ret_t, Eigen::Dynamic, 1> g(3);
  g << 0.0, 0.0, 0.0;

  Eigen::Array<ret_t, Eigen::Dynamic, 1> log_g_old(3);
  log_g_old << NEGATIVE_INFTY, NEGATIVE_INFTY, NEGATIVE_INFTY;

  ret_t log_t_old = 0.0;
  ret_t log_t_new = 0.0;
  int sign_z = sign(z_t);
  auto log_z = log(abs(z_t));

  double log_precision = log(precision);
  int log_t_new_sign = 1.0;
  int log_t_old_sign = 1.0;

  Eigen::ArrayXi log_g_old_sign(3);
  log_g_old_sign << 1, 1, 1;

  int sign_zk = sign_z;
  int k = 0;
  ret_t inner_diff = 1;
  Eigen::Matrix<ret_t, -1, 1> g_current(3);
  g_current << 0, 0, 0;

  while (inner_diff > precision && k < max_steps) {
    ret_t p = ((a1_t + k) * (a2_t + k) / ((b1 + k) * (1.0 + k)));
    if (p == 0) {
      return;
    }
    log_t_new += log(fabs(p)) + log_z;
    log_t_new_sign = sign(value_of_rec(p)) * log_t_new_sign;

    if (calc_a1_euler) {
      ret_t term_a1
          = log_g_old_sign(0) * log_t_old_sign * exp(log_g_old(0) - log_t_old)
            + inv(a1_t + k);
      log_g_old(0) = log_t_new + log(abs(term_a1));
      log_g_old_sign(0) = sign(value_of_rec(term_a1)) * log_t_new_sign;
      g_current(0) = log_g_old_sign(0) * exp(log_g_old(0)) * sign_zk;
      g(0) += g_current(0);
    }

    if (calc_a2_euler) {
      ret_t term_a2
          = log_g_old_sign(1) * log_t_old_sign * exp(log_g_old(1) - log_t_old)
            + inv(a2_t + k);
      log_g_old(1) = log_t_new + log(abs(term_a2));
      log_g_old_sign(1) = sign(value_of_rec(term_a2)) * log_t_new_sign;
      g_current(1) = log_g_old_sign(1) * exp(log_g_old(1)) * sign_zk;
      g(1) += g_current(1);
    }

    if (calc_b1) {
      ret_t term_b1
          = log_g_old_sign(2) * log_t_old_sign * exp(log_g_old(2) - log_t_old)
            + inv(-(b1 + k));
      log_g_old(2) = log_t_new + log(abs(term_b1));
      log_g_old_sign(2) = sign(value_of_rec(term_b1)) * log_t_new_sign;
      g_current(2) = log_g_old_sign(2) * exp(log_g_old(2)) * sign_zk;
      g(2) += g_current(2);
    }

    inner_diff = g_current.array().abs().maxCoeff();

    log_t_old = log_t_new;
    log_t_old_sign = log_t_new_sign;
    sign_zk *= sign_z;
    ++k;
  }

  auto pre_mult_ab = inv(pow(1.0 - z, a2));

  if (calc_a1) {
    g_a1 = euler_transform ? -pre_mult_ab * g(1) : g(0);
  }

  if (calc_a2) {
    if (euler_transform) {
      auto hyper_da2 = hypergeometric_2F1(a1_t, a2, b1, z_t);
      g_a2 = -pre_mult_ab * hyper_da2 * log1m(z) + pre_mult_ab * g(0);
    } else {
      g_a2 = g(1);
    }
  }

  if (calc_b1) {
    g_b1 = euler_transform ? pre_mult_ab * (g(1) + g(2)) : g(2);
  }

  if (k > max_steps) {
    throw_domain_error("grad_2F1", "k (internal counter)", max_steps,
                       "exceeded ",
                       " iterations, hypergeometric function gradient "
                       "did not converge.");
  }
  return;
}
}  // namespace internal

/**
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * @tparam Tg1 type of outcome variable to store a1 gradient
 * @tparam Tg2 type of outcome variable to store a2 gradient
 * @tparam Tg3 type of outcome variable to store b1 gradient
 * @tparam T_gz type of outcome variable to store z gradient
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
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
template <typename Tg1, typename Tg2, typename Tg3, typename T_gz, typename T1,
          typename T2, typename T3, typename T_z>
void grad_2F1(Tg1& g_a1, Tg2& g_a2, Tg3& g_b1, T_gz& g_z, const T1& a1,
              const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  internal::grad_2F1_impl<!is_constant<T1>::value, !is_constant<T2>::value,
                          !is_constant<T3>::value, !is_constant<T_z>::value>(
      g_a1, g_a2, g_b1, g_z, value_of(a1), value_of(a2), value_of(b1),
      value_of(z), precision, max_steps);
  return;
}

/**
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * Overload for use where the destination gradients should be the same type
 * as the input variables (needed for the grad_inc_beta overloads)
 *
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
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
template <typename T1, typename T2, typename T3, typename T_z>
void grad_2F1(T1& g_a1, T2& g_a2, T3& g_b1, T_z& g_z, const T1& a1,
              const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  internal::grad_2F1_impl<true, true, true, true>(g_a1, g_a2, g_b1, g_z, a1, a2,
                                                  b1, z, precision, max_steps);
  return;
}

}  // namespace math
}  // namespace stan
#endif
