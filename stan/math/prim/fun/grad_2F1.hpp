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
 * function, 2F1 with respect to the a1, a2, and b2 parameters.
 *
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * @tparam calc_a1 boolean for whether to calculate gradients w.r.t a1
 * @tparam calc_a2 boolean for whether to calculate gradients w.r.t a2
 * @tparam calc_b1 boolean for whether to calculate gradients w.r.t b1
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 * @return Three-element tuple containing gradients w.r.t. a1, a2, and b1,
 *   as indicated by the calc_a1, calc_a2, and calc_b1 booleans
 */
template <bool calc_a1, bool calc_a2, bool calc_b1, typename T1, typename T2,
          typename T3, typename T_z,
          typename ScalarT = return_type_t<T1, T2, T3, T_z>,
          typename TupleT = std::tuple<ScalarT, ScalarT, ScalarT>>
TupleT grad_2F1_impl_ab(const T1& a1, const T2& a2, const T3& b1, const T_z& z,
                        double precision = 1e-14, int max_steps = 1e6) {
  TupleT grad_tuple = TupleT(0, 0, 0);

  if (z == 0) {
    return grad_tuple;
  }

  using ScalarArrayT = Eigen::Array<ScalarT, 3, 1>;
  ScalarArrayT log_g_old = ScalarArrayT::Constant(3, 1, NEGATIVE_INFTY);

  ScalarT log_t_old = 0.0;
  ScalarT log_t_new = 0.0;
  int sign_z = sign(z);
  auto log_z = log(abs(z));

  double log_precision = log(precision);
  int log_t_new_sign = 1.0;
  int log_t_old_sign = 1.0;

  Eigen::Array<int, 3, 1> log_g_old_sign = Eigen::Array<int, 3, 1>::Ones(3);

  int sign_zk = sign_z;
  int k = 0;
  const int min_steps = 5;
  ScalarT inner_diff = 1;
  ScalarArrayT g_current = ScalarArrayT::Zero(3);

  while ((inner_diff > precision || k < min_steps) && k < max_steps) {
    ScalarT p = ((a1 + k) * (a2 + k) / ((b1 + k) * (1.0 + k)));
    if (p == 0) {
      return grad_tuple;
    }
    log_t_new += log(fabs(p)) + log_z;
    log_t_new_sign = sign(value_of_rec(p)) * log_t_new_sign;

    if (calc_a1) {
      ScalarT term_a1
          = log_g_old_sign(0) * log_t_old_sign * exp(log_g_old(0) - log_t_old)
            + inv(a1 + k);
      log_g_old(0) = log_t_new + log(abs(term_a1));
      log_g_old_sign(0) = sign(value_of_rec(term_a1)) * log_t_new_sign;
      g_current(0) = log_g_old_sign(0) * exp(log_g_old(0)) * sign_zk;
      std::get<0>(grad_tuple) += g_current(0);
    }

    if (calc_a2) {
      ScalarT term_a2
          = log_g_old_sign(1) * log_t_old_sign * exp(log_g_old(1) - log_t_old)
            + inv(a2 + k);
      log_g_old(1) = log_t_new + log(abs(term_a2));
      log_g_old_sign(1) = sign(value_of_rec(term_a2)) * log_t_new_sign;
      g_current(1) = log_g_old_sign(1) * exp(log_g_old(1)) * sign_zk;
      std::get<1>(grad_tuple) += g_current(1);
    }

    if (calc_b1) {
      ScalarT term_b1
          = log_g_old_sign(2) * log_t_old_sign * exp(log_g_old(2) - log_t_old)
            + inv(-(b1 + k));
      log_g_old(2) = log_t_new + log(abs(term_b1));
      log_g_old_sign(2) = sign(value_of_rec(term_b1)) * log_t_new_sign;
      g_current(2) = log_g_old_sign(2) * exp(log_g_old(2)) * sign_zk;
      std::get<2>(grad_tuple) += g_current(2);
    }

    inner_diff = g_current.array().abs().maxCoeff();

    log_t_old = log_t_new;
    log_t_old_sign = log_t_new_sign;
    sign_zk *= sign_z;
    ++k;
  }

  if (k > max_steps) {
    throw_domain_error("grad_2F1", "k (internal counter)", max_steps,
                       "exceeded ",
                       " iterations, hypergeometric function gradient "
                       "did not converge.");
  }
  return grad_tuple;
}

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
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 * @return Four-element tuple containing gradients w.r.t. to each parameter,
 *   as indicated by the calc_* booleans
 */
template <bool calc_a1, bool calc_a2, bool calc_b1, bool calc_z, typename T1,
          typename T2, typename T3, typename T_z,
          typename ScalarT = return_type_t<T1, T2, T3, T_z>,
          typename TupleT = std::tuple<ScalarT, ScalarT, ScalarT, ScalarT>>
TupleT grad_2F1_impl(const T1& a1, const T2& a2, const T3& b1, const T_z& z,
                     double precision = 1e-14, int max_steps = 1e6) {
  bool euler_transform = false;
  try {
    check_2F1_converges("hypergeometric_2F1", a1, a2, b1, z);
  } catch (const std::exception& e) {
    // Apply Euler's hypergeometric transformation if function
    // will not converge with current arguments
    check_2F1_converges("hypergeometric_2F1 (euler transform)", b1 - a1, a2, b1,
                        z / (z - 1));
    euler_transform = true;
  }

  std::tuple<ScalarT, ScalarT, ScalarT> grad_tuple_ab;
  TupleT grad_tuple_rtn = TupleT(0, 0, 0, 0);
  if (euler_transform) {
    ScalarT a1_euler = a2;
    ScalarT a2_euler = b1 - a1;
    ScalarT z_euler = z / (z - 1);
    if (calc_z) {
      auto hyper1 = hypergeometric_2F1(a1_euler, a2_euler, b1, z_euler);
      auto hyper2 = hypergeometric_2F1(1 + a2, 1 - a1 + b1, 1 + b1, z_euler);
      auto pre_mult = a2 * pow(1 - z, -1 - a2);
      std::get<3>(grad_tuple_rtn)
          = a2 * pow(1 - z, -1 - a2) * hyper1
            + (a2 * (b1 - a1) * pow(1 - z, -a2)
               * (inv(z - 1) - z / square(z - 1)) * hyper2)
                  / b1;
    }
    if (calc_a1 || calc_a2 || calc_b1) {
      // 'a' gradients under Euler transform are constructed using the gradients
      // of both elements, so need to compute both if any are required
      constexpr bool calc_a1_euler = calc_a1 || calc_a2;
      // 'b' gradients under Euler transform require gradients from 'a2'
      constexpr bool calc_a2_euler = calc_a1 || calc_a2 || calc_b1;
      grad_tuple_ab = grad_2F1_impl_ab<calc_a1_euler, calc_a2_euler, calc_b1>(
          a1_euler, a2_euler, b1, z_euler);

      auto pre_mult_ab = inv(pow(1.0 - z, a2));
      if (calc_a1) {
        std::get<0>(grad_tuple_rtn) = -pre_mult_ab * std::get<1>(grad_tuple_ab);
      }
      if (calc_a2) {
        auto hyper_da2 = hypergeometric_2F1(a1_euler, a2, b1, z_euler);
        std::get<1>(grad_tuple_rtn)
            = -pre_mult_ab * hyper_da2 * log1m(z)
              + pre_mult_ab * std::get<0>(grad_tuple_ab);
      }
      if (calc_b1) {
        std::get<2>(grad_tuple_rtn)
            = pre_mult_ab
              * (std::get<1>(grad_tuple_ab) + std::get<2>(grad_tuple_ab));
      }
    }
  } else {
    if (calc_z) {
      auto hyper_2f1_dz = hypergeometric_2F1(a1 + 1.0, a2 + 1.0, b1 + 1.0, z);
      std::get<3>(grad_tuple_rtn) = (a1 * a2 * hyper_2f1_dz) / b1;
    }
    if (calc_a1 || calc_a2 || calc_b1) {
      grad_tuple_ab
          = grad_2F1_impl_ab<calc_a1, calc_a2, calc_b1>(a1, a2, b1, z);
      if (calc_a1) {
        std::get<0>(grad_tuple_rtn) = std::get<0>(grad_tuple_ab);
      }
      if (calc_a2) {
        std::get<1>(grad_tuple_rtn) = std::get<1>(grad_tuple_ab);
      }
      if (calc_b1) {
        std::get<2>(grad_tuple_rtn) = std::get<2>(grad_tuple_ab);
      }
    }
  }
  return grad_tuple_rtn;
}
}  // namespace internal

/**
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * Overload for use where the destination gradients are not required to be the
 * same type as the input variables (most use-cases except grad_inc_beta)
 *
 * @tparam ReturnSameT Whether return gradients need to be the same type as
 *  as inputs
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 */
template <bool ReturnSameT, typename T1, typename T2, typename T3, typename T_z,
          require_not_t<std::integral_constant<bool, ReturnSameT>>* = nullptr>
auto grad_2F1(const T1& a1, const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  return internal::grad_2F1_impl<
      !is_constant<T1>::value, !is_constant<T2>::value, !is_constant<T3>::value,
      !is_constant<T_z>::value>(value_of(a1), value_of(a2), value_of(b1),
                                value_of(z), precision, max_steps);
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
 * @tparam ReturnSameT Whether return gradients need to be the same type as
 *  as inputs
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 */
template <bool ReturnSameT, typename T1, typename T2, typename T3, typename T_z,
          require_t<std::integral_constant<bool, ReturnSameT>>* = nullptr>
auto grad_2F1(const T1& a1, const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  return internal::grad_2F1_impl<true, true, true, true>(a1, a2, b1, z,
                                                         precision, max_steps);
}

/**
 * Calculate the gradients of the hypergeometric function (2F1)
 * as the power series stopping when the series converges
 * to within <code>precision</code> or throwing when the
 * function takes <code>max_steps</code> steps.
 *
 * @tparam T1 scalar type of a1
 * @tparam T2 scalar type of a2
 * @tparam T3 scalar type of b1
 * @tparam T_z scalar type of z
 * @param[in] a1 see generalized hypergeometric function definition
 * @param[in] a2 see generalized hypergeometric function definition
 * @param[in] b1 see generalized hypergeometric function definition
 * @param[in] z see generalized hypergeometric function definition
 * @param[in] precision magnitude of the increment of the infinite sum
 *   to truncate the sum
 * @param[in] max_steps number of steps to take
 */
template <typename T1, typename T2, typename T3, typename T_z>
auto grad_2F1(const T1& a1, const T2& a2, const T3& b1, const T_z& z,
              double precision = 1e-14, int max_steps = 1e6) {
  return grad_2F1<false>(a1, a2, b1, z, precision, max_steps);
}

}  // namespace math
}  // namespace stan
#endif
