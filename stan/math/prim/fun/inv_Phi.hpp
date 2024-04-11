#ifndef STAN_MATH_PRIM_FUN_INV_PHI_HPP
#define STAN_MATH_PRIM_FUN_INV_PHI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
/**
 *  The largest integer that protects against floating point errors
 * for the inv_Phi function. The value was found by finding the largest
 * integer that passed the unit tests for accuracy when the input into inv_Phi
 * is near 1.
 */
const int BIGINT = 2000000000;

/**
 * The inverse of the unit normal cumulative distribution function.
 *
 * @tparam LogP Whether the input probability is already on the log scale.
 * @param p argument between 0 and 1 inclusive
 * @return Real value of the inverse cdf for the standard normal distribution.
 */
template <bool LogP = false>
inline double inv_Phi_impl(double p) {
  static constexpr double log_a[8]
      = {1.2199838032983212, 4.8914137334471356, 7.5865960847956080,
         9.5274618535358388, 10.734698580862359, 11.116406781896242,
         10.417226196842595, 7.8276718012189362};
  static constexpr double log_b[8] = {0.,
                                      3.7451021830139207,
                                      6.5326064640478618,
                                      8.5930788436817044,
                                      9.9624069236663077,
                                      10.579180688621286,
                                      10.265665328832871,
                                      8.5614962136628454};
  static constexpr double log_c[8]
      = {0.3530744474482423, 1.5326298343683388, 1.7525849400614634,
         1.2941374937060454, 0.2393776640901312, -1.419724057885092,
         -3.784340465764968, -7.163234779359426};
  static constexpr double log_d[8] = {0.0,
                                      0.71939547349472054982,
                                      0.51663958798453168964,
                                      -0.37140093392784434556,
                                      -1.9098407084572139869,
                                      -4.186547581055928724,
                                      -7.5099767712254150709,
                                      -20.673761573859248841};
  static constexpr double log_e[8]
      = {1.8958048169567149, 1.6981417567726154, 0.5793212339927351,
         -1.215503791936417, -3.629396584023968, -6.690500273261249,
         -10.51540298415323, -15.41979457491781};
  static constexpr double log_f[8] = {0.,
                                      -0.511105318617135,
                                      -1.988286302259815,
                                      -4.208049039384857,
                                      -7.147448611626374,
                                      -10.89973190740069,
                                      -15.76637472711685,
                                      -33.82373901099482};

  double log_p = LogP ? p : log(p);

  double log_q = log_p <= LOG_HALF ? log_diff_exp(LOG_HALF, log_p)
                                   : log_diff_exp(log_p, LOG_HALF);
  int log_q_sign = log_p <= LOG_HALF ? -1 : 1;
  double log_r = log_q_sign == -1 ? log_p : log1m_exp(log_p);

  if (stan::math::is_inf(log_r)) {
    return 0;
  }

  double log_inner_r;
  double log_pre_mult;
  const double* num_ptr;
  const double* den_ptr;

  static constexpr double LOG_FIVE = LOG_TEN - LOG_TWO;
  static constexpr double LOG_16 = LOG_TWO * 4;
  static constexpr double LOG_425 = 6.0520891689244171729;
  static constexpr double LOG_425_OVER_1000 = LOG_425 - LOG_TEN * 3;

  if (log_q <= LOG_425_OVER_1000) {
    log_inner_r = log_diff_exp(LOG_425_OVER_1000 * 2, log_q * 2);
    log_pre_mult = log_q;
    num_ptr = &log_a[0];
    den_ptr = &log_b[0];
  } else {
    double log_temp_r = log(-log_r) / 2.0;
    if (log_temp_r <= LOG_FIVE) {
      log_inner_r = log_diff_exp(log_temp_r, LOG_16 - LOG_TEN);
      num_ptr = &log_c[0];
      den_ptr = &log_d[0];
    } else {
      log_inner_r = log_diff_exp(log_temp_r, LOG_FIVE);
      num_ptr = &log_e[0];
      den_ptr = &log_f[0];
    }
    log_pre_mult = 0.0;
  }

  // As computation requires evaluating r^8, this causes a loss of precision,
  // even when on the log space. We can mitigate this by scaling the
  // exponentiated result (dividing by 10), since the same scaling is applied
  // to the numerator and denominator.
  Eigen::VectorXd log_r_pow
      = Eigen::ArrayXd::LinSpaced(8, 0, 7) * log_inner_r - LOG_TEN;
  Eigen::Map<const Eigen::VectorXd> num_map(num_ptr, 8);
  Eigen::Map<const Eigen::VectorXd> den_map(den_ptr, 8);
  double log_result
      = log_sum_exp(log_r_pow + num_map) - log_sum_exp(log_r_pow + den_map);
  return log_q_sign * exp(log_pre_mult + log_result);
}
}  // namespace internal

/**
 * Return the value of the inverse standard normal cumulative distribution
 * function at the specified argument.
 *
 * The precision is at or better than 1.5e-15 for values between 0.0000001 he
 * largest integer that protects against floating point errors for the inv_Phi
 * function. The value was found by finding the largest integer that passed the
 * unit tests for accuracy when the input into inv_Phi is near 1.
 *
 * @param p argument between 0 and 1 inclusive
 * @return real value of the inverse cdf for the standard normal distribution
 */
inline double inv_Phi(double p) {
  check_bounded("inv_Phi", "Probability variable", p, 0, 1);

  if (p < 8e-311) {
    return NEGATIVE_INFTY;
  }
  if (p == 1) {
    return INFTY;
  }
  return p >= 0.9999 ? -internal::inv_Phi_impl(
             (internal::BIGINT - internal::BIGINT * p) / internal::BIGINT)
                     : internal::inv_Phi_impl(p);
}

/**
 * Structure to wrap inv_Phi() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable in range [0, 1]
 * @return Inverse unit normal CDF of x.
 * @throw std::domain_error if x is not between 0 and 1.
 */
struct inv_Phi_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return inv_Phi(x);
  }
};

/**
 * Vectorized version of inv_Phi().
 *
 * @tparam T type of container
 * @param x variables in range [0, 1]
 * @return Inverse unit normal CDF of each value in x.
 * @throw std::domain_error if any value is not between 0 and 1.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto inv_Phi(const T& x) {
  return apply_scalar_unary<inv_Phi_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
