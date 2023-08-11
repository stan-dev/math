#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LOG_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse of the unit normal cumulative distribution function evaluated at
 * the log probability.
 *
 * @param log_p argument between -Inf and 0 inclusive
 * @return Real value of the inverse cdf for the standard normal distribution.
 */
inline double std_normal_log_qf(double log_p) {
  check_not_nan("std_normal_log_qf", "Log probability variable", log_p);
  check_less_or_equal("std_normal_log_qf", "Probability variable", log_p, 0);

  if (log_p == NEGATIVE_INFTY) {
    return NEGATIVE_INFTY;
  }
  if (log_p == 0) {
    return INFTY;
  }

  static const double log_a[8]
      = {1.2199838032983212, 4.8914137334471356, 7.5865960847956080,
         9.5274618535358388, 10.734698580862359, 11.116406781896242,
         10.417226196842595, 7.8276718012189362};
  static const double log_b[8] = {0.,
                                  3.7451021830139207,
                                  6.5326064640478618,
                                  8.5930788436817044,
                                  9.9624069236663077,
                                  10.579180688621286,
                                  10.265665328832871,
                                  8.5614962136628454};
  static const double log_c[8]
      = {0.3530744474482423, 1.5326298343683388, 1.7525849400614634,
         1.2941374937060454, 0.2393776640901312, -1.419724057885092,
         -3.784340465764968, -7.163234779359426};
  static const double log_d[8] = {0.,
                                  0.7193954734947205,
                                  0.5166395879845317,
                                  -0.371400933927844,
                                  -1.909840708457214,
                                  -4.186547581055928,
                                  -7.509976771225415,
                                  -20.67376157385924};
  static const double log_e[8]
      = {1.8958048169567149, 1.6981417567726154, 0.5793212339927351,
         -1.215503791936417, -3.629396584023968, -6.690500273261249,
         -10.51540298415323, -15.41979457491781};
  static const double log_f[8] = {0.,
                                  -0.511105318617135,
                                  -1.988286302259815,
                                  -4.208049039384857,
                                  -7.147448611626374,
                                  -10.89973190740069,
                                  -15.76637472711685,
                                  -33.82373901099482};

  double val;
  double log_q = log_p <= LOG_HALF ? log_diff_exp(LOG_HALF, log_p)
                                   : log_diff_exp(log_p, LOG_HALF);
  int log_q_sign = log_p <= LOG_HALF ? -1 : 1;

  if (log_q <= -0.85566611005772) {
    double log_r = log_diff_exp(-1.71133222011544, 2 * log_q);
    double log_agg_a = log_sum_exp(log_a[7] + log_r, log_a[6]);
    double log_agg_b = log_sum_exp(log_b[7] + log_r, log_b[6]);

    for (int i = 0; i < 6; i++) {
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[5 - i]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[5 - i]);
    }

    return log_q_sign * exp(log_q + log_agg_a - log_agg_b);
  } else {
    double log_r = log_q_sign == -1 ? log_p : log1m_exp(log_p);

    if (stan::math::is_inf(log_r)) {
      return 0;
    }

    log_r = log(sqrt(-log_r));

    if (log_r <= 1.60943791243410) {
      log_r = log_diff_exp(log_r, 0.47000362924573);
      double log_agg_c = log_sum_exp(log_c[7] + log_r, log_c[6]);
      double log_agg_d = log_sum_exp(log_d[7] + log_r, log_d[6]);

      for (int i = 0; i < 6; i++) {
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[5 - i]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[5 - i]);
      }
      val = exp(log_agg_c - log_agg_d);
    } else {
      log_r = log_diff_exp(log_r, 1.60943791243410);
      double log_agg_e = log_sum_exp(log_e[7] + log_r, log_e[6]);
      double log_agg_f = log_sum_exp(log_f[7] + log_r, log_f[6]);

      for (int i = 0; i < 6; i++) {
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[5 - i]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[5 - i]);
      }
      val = exp(log_agg_e - log_agg_f);
    }
    if (log_q_sign == -1)
      return -val;
  }
  return val;
}

/**
 * Structure to wrap std_normal_log_qf() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable in range [-Inf, 0]
 * @return inverse of the unit normal CDF of x
 * @throw std::domain_error if x is not less than or equal to 0
 */
struct std_normal_log_qf_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return std_normal_log_qf(x);
  }
};

/**
 * A vectorized version of std_normal_log_qf() that accepts std::vectors, Eigen
 * Matrix/Array objects, or expressions, and containers of these.
 *
 * @tparam T type of container
 * @param x container
 * @return inverse unit normal CDF of each value in x
 * @throw std::domain_error if x is not less than or equal to 0
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto std_normal_log_qf(const T& x) {
  return apply_scalar_unary<std_normal_log_qf_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
